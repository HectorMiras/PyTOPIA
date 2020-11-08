function [] = calibracion(D)

%% CORREGIR INPUT DE LAS DOSIS
D(D==0)=[]; %Quitar si se ha introducido el valor de dosis 0
D=reshape(D,[length(D),1]); %El valor de dosis tiene que tener formato de columna

N_images=1;
N_data_image=length(D);

%% SELECCIONAR DATOS DE LAS PELÍCULAS

xmin0=zeros(1,N_images);
xmax0=zeros(1,N_images);
ymin0=zeros(1,N_images);
ymax0=zeros(1,N_images);

xmin=zeros(1,length(D));
xmax=zeros(1,length(D));
ymin=zeros(1,length(D));
ymax=zeros(1,length(D));

dat_archivo=cell(1,N_images);

count=1;

for zz=1:N_images

% Buscamos por explorador el archivo
[archivo,ruta]=uigetfile('*.tif','Abrir un archivo de datos');

% Pegamos el nombre del archivo y de la dirección
if archivo==0
 return;
else
 dat_archivo(1,zz)=cellstr(strcat(ruta,archivo));
end

I = imread(char(dat_archivo(1,zz)));

% SIN IRRADIAR

% Sección de la imagen donde se extrae la información
h=figure;
imshow(I), axis('image'); title('SELECCIONE ÁREA DE DOSIS D=0Gy');
sp=getrect(h);
close(h)

% Get the x and y co-ordinates
xmin0(zz) = max(floor(sp(1)), 1);
ymin0(zz) = max(floor(sp(2)), 1);
xmax0(zz) = min(ceil(sp(1) + sp(3)));
ymax0(zz) = min(ceil(sp(2) +sp(4)));

% IRRADIADO
for ii=count:count-1+N_data_image(zz)
% Sección de la imagen donde se extrae la información
valor_D=num2str(D(ii));
h=figure;
imshow(I), axis('image'); title(strcat('SELECCIONE ÁREA DE DOSIS D=',valor_D,'Gy'));
sp=getrect(h);
close(h)

% Get the x and y co-ordinates
xmin(ii) = max(floor(sp(1)), 1);
ymin(ii) = max(floor(sp(2)), 1);
xmax(ii) = min(ceil(sp(1) + sp(3)));
ymax(ii) = min(ceil(sp(2) +sp(4)));

end

count=count+N_data_image(zz);

end

N_x=min([xmax0-xmin0+1, xmax-xmin+1]); % Número de columnas del cuadrado
N_y=min([ymax0-ymin0+1, ymax-ymin+1]); % Número de filas del cuadrado
N_data=N_x*N_y; % NOTA: Se le suma +1 porque hay que contar el primer pixel

%% OBTENCIÓN DE LOS DATOS DE LA PELÍCULA

I_R=zeros(length(D),N_data);
I_G=zeros(length(D),N_data);
I_B=zeros(length(D),N_data);

I_R0=zeros(N_images,N_data);
I_G0=zeros(N_images,N_data);
I_B0=zeros(N_images,N_data);

count=1;

for zz=1:N_images

I = imread(char(dat_archivo(1,zz)));

% SIN IRRADIAR
I_R0(zz,:)=reshape(I(ymin0(zz):ymin0(zz)+N_y-1,xmin0(zz):xmin0(zz)+N_x-1,1),1,[]);
I_G0(zz,:)=reshape(I(ymin0(zz):ymin0(zz)+N_y-1,xmin0(zz):xmin0(zz)+N_x-1,2),1,[]);
I_B0(zz,:)=reshape(I(ymin0(zz):ymin0(zz)+N_y-1,xmin0(zz):xmin0(zz)+N_x-1,3),1,[]);

% IRRADIADO
for ii=count:count-1+N_data_image(zz)
I_R(ii,:)=reshape(I(ymin(ii):ymin(ii)+N_y-1,xmin(ii):xmin(ii)+N_x-1,1),1,[]);
I_G(ii,:)=reshape(I(ymin(ii):ymin(ii)+N_y-1,xmin(ii):xmin(ii)+N_x-1,2),1,[]);
I_B(ii,:)=reshape(I(ymin(ii):ymin(ii)+N_y-1,xmin(ii):xmin(ii)+N_x-1,3),1,[]);
end

count=count+N_data_image(zz);

end

clearvars -except D I_R I_G I_B I_R0 I_G0 I_B0 N_images N_data_image N_x N_y N_data

%% CÁLCULOS DATOS NOD PARA CURVA DE CALIBRACIÓN

% MEDIA DE LOS PÍXELES Y netOD
I_R0_media=sum(I_R0,2)/N_data;
I_G0_media=sum(I_G0,2)/N_data;
I_B0_media=sum(I_B0,2)/N_data;

I_R_media=sum(I_R,2)/N_data;
I_G_media=sum(I_G,2)/N_data;
I_B_media=sum(I_B,2)/N_data;

sigma_I_R0_media=sqrt(sum((I_R0-I_R0_media).^2,2)/(N_data-1));
sigma_I_G0_media=sqrt(sum((I_G0-I_G0_media).^2,2)/(N_data-1));
sigma_I_B0_media=sqrt(sum((I_B0-I_B0_media).^2,2)/(N_data-1));

sigma_I_R_media=sqrt(sum((I_R-I_R_media).^2,2)/(N_data-1));
sigma_I_G_media=sqrt(sum((I_G-I_G_media).^2,2)/(N_data-1));
sigma_I_B_media=sqrt(sum((I_B-I_B_media).^2,2)/(N_data-1));

count=1;
for zz=1:N_images
I_R_media_netOD(count:count+N_data_image(zz)-1,:)=-log10(I_R_media(count:count+N_data_image(zz)-1,:)./I_R0_media(zz));
I_G_media_netOD(count:count+N_data_image(zz)-1,:)=-log10(I_G_media(count:count+N_data_image(zz)-1,:)./I_G0_media(zz));
I_B_media_netOD(count:count+N_data_image(zz)-1,:)=-log10(I_B_media(count:count+N_data_image(zz)-1,:)./I_B0_media(zz));

sigma_I_R_media_netOD(count:count+N_data_image(zz)-1,:)=sqrt((sigma_I_R_media(count:count+N_data_image(zz)-1,:)./I_R_media(count:count+N_data_image(zz)-1,:)).^2+(sigma_I_R0_media(zz)./I_R0_media(zz)).^2)/log(10); 
sigma_I_G_media_netOD(count:count+N_data_image(zz)-1,:)=sqrt((sigma_I_G_media(count:count+N_data_image(zz)-1,:)./I_G_media(count:count+N_data_image(zz)-1,:)).^2+(sigma_I_G0_media(zz)./I_G0_media(zz)).^2)/log(10);
sigma_I_B_media_netOD(count:count+N_data_image(zz)-1,:)=sqrt((sigma_I_B_media(count:count+N_data_image(zz)-1,:)./I_B_media(count:count+N_data_image(zz)-1,:)).^2+(sigma_I_B0_media(zz)./I_B0_media(zz)).^2)/log(10);

count=count+N_data_image(zz);
end

% netOD Y MEDIA
count=1;
for zz=1:N_images
netOD_R(count:count+N_data_image(zz)-1,:)=-log10(I_R(count:count+N_data_image(zz)-1,:)./I_R0_media(zz));
netOD_G(count:count+N_data_image(zz)-1,:)=-log10(I_G(count:count+N_data_image(zz)-1,:)./I_G0_media(zz));
netOD_B(count:count+N_data_image(zz)-1,:)=-log10(I_B(count:count+N_data_image(zz)-1,:)./I_B0_media(zz));

count=count+N_data_image(zz);
end

netOD_R_media=sum(netOD_R,2)/N_data;
netOD_G_media=sum(netOD_G,2)/N_data;
netOD_B_media=sum(netOD_B,2)/N_data;

sigma_netOD_R_media=sqrt(sum((netOD_R-netOD_R_media).^2,2)/(N_data-1));
sigma_netOD_G_media=sqrt(sum((netOD_G-netOD_G_media).^2,2)/(N_data-1));
sigma_netOD_B_media=sqrt(sum((netOD_B-netOD_B_media).^2,2)/(N_data-1));

% PESOS
weights_R_NODm=1./(sigma_I_R_media_netOD).^2; weights_R_NODm=weights_R_NODm./sum(weights_R_NODm); %Usando el NOD de los valores medios
weights_G_NODm=1./(sigma_I_G_media_netOD).^2; weights_G_NODm=weights_G_NODm./sum(weights_G_NODm);
weights_B_NODm=1./(sigma_I_B_media_netOD).^2; weights_B_NODm=weights_B_NODm./sum(weights_B_NODm);

weights_R_mNOD=1./(sigma_netOD_R_media).^2; weights_R_mNOD=weights_R_mNOD./sum(weights_R_mNOD); %Usando la media de la NOD
weights_G_mNOD=1./(sigma_netOD_G_media).^2; weights_G_mNOD=weights_G_mNOD./sum(weights_G_mNOD);
weights_B_mNOD=1./(sigma_netOD_B_media).^2; weights_B_mNOD=weights_B_mNOD./sum(weights_B_mNOD);

%% PRECALIBRACIÓN (hallar y fijar el valor óptimo de n)

% OPCIONES GENERALES DEL AJUSTE
options = fitoptions('Method','NonlinearLeastSquares');
options.Algorithm = 'Trust-Region';
options.StartPoint = [ 5 50 2];
options.lower = [0 0 0];
options.DiffMaxChange = 0.1;
options.DiffMinChange = 1e-8;
options.MaxFunEvals = 600;
options.MaxIter = 500;
options.TolFun = 1e-6;
options.TolX = 1e-6;

% FUNCIÓN DE AJUSTE
calibration=@(b,c,n,x) b.*x+c.*x.^n;
    
% AJUSTE: NOD de los valores medios

% Red
options.Weights = weights_R_NODm;
[cal_object,~] = fit( I_R_media_netOD, D, calibration, options);
coef_R_NODm = coeffvalues(cal_object); coef_R_NODm(3)=round(coef_R_NODm(3)*100)/100;
% Green
options.Weights = weights_G_NODm;
[cal_object,~] = fit( I_G_media_netOD, D, calibration, options);
coef_G_NODm = coeffvalues(cal_object); coef_G_NODm(3)=round(coef_G_NODm(3)*100)/100;

% Blue
options.Weights = weights_B_NODm;
[cal_object,~] = fit( I_B_media_netOD, D, calibration, options);
coef_B_NODm = coeffvalues(cal_object); coef_B_NODm(3)=round(coef_B_NODm(3)*100)/100;

% AJUSTE: valores medios de los NOD

% Red
options.Weights = weights_R_mNOD;
[cal_object,~] = fit( netOD_R_media, D, calibration, options);
coef_R_mNOD = coeffvalues(cal_object); coef_R_mNOD(3) = round(coef_R_mNOD(3)*100)/100;

% Green
options.Weights = weights_G_mNOD;
[cal_object,~] = fit( netOD_G_media, D, calibration, options);
coef_G_mNOD = coeffvalues(cal_object); coef_G_mNOD(3) = round(coef_G_mNOD(3)*100)/100;

% Blue
options.Weights = weights_B_mNOD;
[cal_object,~] = fit( netOD_B_media, D, calibration, options);
coef_B_mNOD = coeffvalues(cal_object); coef_B_mNOD(3) = round(coef_B_mNOD(3)*100)/100;


%% CALIBRACIÓN

% OPCIONES GENERALES DEL AJUSTE
options = fitoptions('Method','NonlinearLeastSquares');
options.Algorithm = 'Trust-Region';
options.StartPoint = [ 1 1 ];
options.lower = [ 0 0 ];
options.DiffMaxChange = 0.1;
options.DiffMinChange = 1e-8;
options.MaxFunEvals = 600;
options.MaxIter = 500;
options.TolFun = 1e-6;
options.TolX = 1e-6;

% FUNCIÓN DE AJUSTE
calibration=@(b,c,n,x) b.*x+c.*x.^n;
    
% AJUSTE: NOD de los valores medios

% Red
options.Weights = weights_R_NODm;
[cal_object,~] = fit( I_R_media_netOD, D, @(b,c,x) calibration(b,c,coef_R_NODm(3),x), options);
coef_R_NODm(1:2) = coeffvalues(cal_object);

% Green
options.Weights = weights_G_NODm;
[cal_object,~] = fit( I_G_media_netOD, D, @(b,c,x) calibration(b,c,coef_G_NODm(3),x), options);
coef_G_NODm(1:2) = coeffvalues(cal_object);

% Blue
options.Weights = weights_B_NODm;
[cal_object,~] = fit( I_B_media_netOD, D, @(b,c,x) calibration(b,c,coef_B_NODm(3),x), options);
coef_B_NODm(1:2) = coeffvalues(cal_object);

% AJUSTE: valores medios de los NOD

% Red
options.Weights = weights_R_mNOD;
[cal_object,~] = fit( netOD_R_media, D, @(b,c,x) calibration(b,c,coef_R_mNOD(3),x), options);
coef_R_mNOD(1:2) = coeffvalues(cal_object);
interval = confint(cal_object,0.95); std_error_R_mNOD=diff(interval)/(2*1.96);

% Green
options.Weights = weights_G_mNOD;
[cal_object,~] = fit( netOD_G_media, D, @(b,c,x) calibration(b,c,coef_G_mNOD(3),x), options);
coef_G_mNOD(1:2) = coeffvalues(cal_object);
interval = confint(cal_object,0.95); std_error_G_mNOD=diff(interval)/(2*1.96);

% Blue
options.Weights = weights_B_mNOD;
[cal_object,~] = fit( netOD_B_media, D, @(b,c,x) calibration(b,c,coef_B_mNOD(3),x), options);
coef_B_mNOD(1:2) = coeffvalues(cal_object);
interval = confint(cal_object,0.95); std_error_B_mNOD=diff(interval)/(2*1.96);

% RESULTADOS
hold on
plot(D,I_R_media_netOD,'rx')
plot(D,I_G_media_netOD,'gx')
plot(D,I_B_media_netOD,'bx')
errorbar(D,netOD_R_media, sigma_netOD_R_media,'rx')
errorbar(D,netOD_G_media, sigma_netOD_G_media,'gx')
errorbar(D,netOD_B_media, sigma_netOD_B_media,'bx')
xlim([0,max(D)])

x=0:0.0001:0.7;
plot(coef_R_NODm(1).*x+coef_R_NODm(2).*x.^coef_R_NODm(3),x,'r-')
plot(coef_G_NODm(1).*x+coef_G_NODm(2).*x.^coef_G_NODm(3),x,'g-')
plot(coef_B_NODm(1).*x+coef_B_NODm(2).*x.^coef_B_NODm(3),x,'b-')
plot(coef_R_mNOD(1).*x+coef_R_mNOD(2).*x.^coef_R_mNOD(3),x,'r:')
plot(coef_G_mNOD(1).*x+coef_G_mNOD(2).*x.^coef_G_mNOD(3),x,'g:')
plot(coef_B_mNOD(1).*x+coef_B_mNOD(2).*x.^coef_B_mNOD(3),x,'b:')

%Renombro y me quedo solo con el ajuste correspondiente al cálculo primero
%del NOD y luego la media.
coef(:,1)=coef_R_mNOD';
coef(:,2)=coef_G_mNOD';
coef(:,3)=coef_B_mNOD';

sigma_coef(:,1)=std_error_R_mNOD';
sigma_coef(:,2)=std_error_G_mNOD';
sigma_coef(:,3)=std_error_B_mNOD';

%% CÁLCULO DE LA MEDIA Y LA DESVIACIÓN DE LOS PARÁMETROS PARA LAS MANCHAS DE CALIBRACIÓN

OD0_cal=mean(log10((2^16-1)./[I_R0_media,I_G0_media,I_B0_media]),1);

netOD_cal=zeros(N_data,3,length(D));
for zz=1:length(D)
    netOD_cal(:,:,zz)=[netOD_R(zz,:); netOD_G(zz,:); netOD_B(zz,:),]';
end

alpha_cal=zeros(N_data,length(D));
beta_cal=zeros(N_data,length(D));

for zz=1:length(D)
    
    Dcal_pixel=coef(1,:).*netOD_cal(:,:,zz)+sign(netOD_cal(:,:,zz)).*coef(2,:).*abs(netOD_cal(:,:,zz)).^(coef(3,:));
    sigma_Dcal=sqrt(netOD_cal(:,:,zz).^2.*sigma_coef(1,:).^2+(netOD_cal(:,:,zz).^2).^(coef(3,:)).*sigma_coef(2,:));

    % Calculo la NOD que me daría la curva de calibración
    NOD=zeros(1,3); % NOD media que debe dar la mancha para la dosis D conocida a partir de la curva de calibración
    for kk=1:3
        NOD(1,kk)=fzero(@(x) coef(1,kk).*x+coef(2,kk).*x.^(coef(3,kk))-D(zz),mean(netOD_cal(:,kk,zz),1));
    end

    d_dose=coef(1,:)+coef(3,:).*coef(2,:).*NOD.^(coef(3,:)-1);
    d_mu_da=d_dose.*NOD;
    d_mu_db=d_dose.*OD0_cal;

    Cai=sum((D(zz)-Dcal_pixel).*d_mu_da./sigma_Dcal.^2,2);
    Caa=sum(d_dose.*NOD.*d_mu_da./sigma_Dcal.^2,2);
    Cab=sum(d_dose.*OD0_cal.*d_mu_da./sigma_Dcal.^2,2);

    Cbi=sum((D(zz)-Dcal_pixel).*d_mu_db./sigma_Dcal.^2,2);
    Cbb=sum(d_dose.*OD0_cal.*d_mu_db./sigma_Dcal.^2,2);

    % Resultado
    alpha_cal(:,zz)=(Cai.*Cbb-Cbi.*Cab)./(Caa.*Cbb-Cab.^2);
    beta_cal(:,zz)=(Cai.*Cab-Cbi.*Caa)./(Cab.^2-Caa.*Cbb);

end

D_cal=[0,reshape(D,[1,length(D)])];

%% GUARDADO DE LOS DATOS
[archivo,ruta]=uiputfile('calibracion.mat','Guardar mapa de dosis sin corregir');

% Pegamos el nombre del archivo y de la dirección
if archivo==0
 return;
else
 dat_archivo=strcat(ruta,archivo);
end

save(dat_archivo, 'coef', 'sigma_coef', 'alpha_cal', 'beta_cal', 'D_cal', 'netOD_cal', 'OD0_cal')
              
end

