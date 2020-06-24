function [] = multicanal

%% CARGA DE LOS COEFICIENTES DE CALIBRACIÓN
[archivo,ruta]=uigetfile('*.mat','Abrir archivo de calibración');

% Pegamos el nombre del archivo y de la dirección
if archivo==0
 return;
else
 dat_archivo=cellstr(strcat(ruta,archivo));
end

% Cargamos el archivo
load(char(dat_archivo));

% Quitar el valor dosis nulo del vector dosis
D_cal(D_cal==0)=[];

%% SELECCIÓN DE LA IMAGEN

% Buscamos por explorador el archivo
[archivo,ruta]=uigetfile('*.tif','Abrir imagen a analizar');

% Pegamos el nombre del archivo y de la dirección
if archivo==0
 return;
else
 dat_archivo=cellstr(strcat(ruta,archivo));
end

I = imread(char(dat_archivo));

[Ny_I,Nx_I,~]=size(I);

%% SELECCIÓN DEL ÁREA SIN IRRADIAR (D=0)

% Sección de la imagen donde se extrae la información
h=figure;
imshow(I), axis('image'), title('SELECIÓN DEL ÁREA NO IRRADIADA');
sp=getrect(h);
close(h)

% Get the x and y co-ordinates
xmin0 = ceil(sp(1)); xmin0(xmin0<=0)=1;
ymin0 = ceil(sp(2)); ymin0(ymin0<=0)=1;
xmax0 = floor(sp(1) + sp(3)); xmax0(xmax0>Nx_I)=Nx_I;
ymax0 = floor(sp(2) + sp(4)); ymax0(ymax0>Ny_I)=Ny_I;

I0=reshape(double(I(ymin0:ymax0,xmin0:xmax0,:)),[(xmax0-xmin0+1)*(ymax0-ymin0+1),3]);

I0_mean=mean(I0,1);

OD0=log10((2^16-1)./I0_mean);

%% SELECCIÓN DEL ÁREA IRRADIADA DE DOSIS CONOCIDA

% Input dosis de referencia
D_ref=input('Introduce valor de dosis de referencia en Gy: ');

% Sección de la imagen donde se extrae la información
h=figure;
imshow(I), axis('image'); title('SELECIÓN DEL ÁREA IRRADIADA DE REFERENCIA');
sp=getrect(h);
close(h)

% Get the x and y co-ordinates
xminD = ceil(sp(1)); xminD(xminD<=0)=1;
yminD = ceil(sp(2)); yminD(yminD<=0)=1;
xmaxD = floor(sp(1) + sp(3)); xmaxD(xmaxD>Nx_I)=Nx_I;
ymaxD = floor(sp(2) + sp(4)); ymaxD(ymaxD>Ny_I)=Ny_I;

ID=reshape(double(I(yminD:ymaxD,xminD:xmaxD,:)),[(xmaxD-xminD+1)*(ymaxD-yminD+1),3]);

netODD=-log10(ID./I0_mean);

netODD_mean=mean(netODD,1);

% Calculo la NOD que me daría la curva de calibración
NOD_D=zeros(1,3); % NOD media que debe dar la mancha para la dosis D conocida a partir de la curva de calibración
for kk=1:3
    NOD_D(1,kk)=fzero(@(x) coef(1,kk).*x+coef(2,kk).*x.^(coef(3,kk))-D_ref,mean(netODD(:,kk),1));
end

% Factor de corrección respecto a la calibración
f_cor=NOD_D./netODD_mean;

% Corrijo los valores a dosis conocida
netODD=f_cor.*netODD;

DD_pixel=(coef(1,:).*netODD+coef(2,:).*netODD.^(coef(3,:)));
sigma_DD=sqrt(netODD.^2.*sigma_coef(1,:).^2+netODD.^(2*coef(3,:)).*sigma_coef(2,:));

%% SELECCIÓN DEL ÁREA IRRADIADA DESCONOCIDA

% Máximo valor de dosis del área desconocida
D_max=input('Introduce el valor máximo de dosis posible en Gy: ');

% Sección de la imagen donde se extrae la información
h=figure;
imshow(I), axis('image'); title('SELECIÓN DEL ÁREA IRRADIADA A ANALIZAR');
sp=getrect(h);
close(h)

% Get the x and y co-ordinates
xmin = ceil(sp(1)); xmin(xmin<=0)=1;
ymin = ceil(sp(2)); ymin(ymin<=0)=1;
xmax = floor(sp(1) + sp(3)); xmax(xmax>Nx_I)=Nx_I;
ymax = floor(sp(2) + sp(4)); ymax(ymax>Ny_I)=Ny_I;
Nx=xmax-xmin+1; Ny=ymax-ymin+1;
N_data = Nx*Ny;
I=reshape(double(I(ymin:ymax,xmin:xmax,:)),[N_data,3]);

netOD=-log10(I./I0_mean); netOD=f_cor.*netOD;

D_pixel=coef(1,:).*netOD+sign(netOD).*coef(2,:).*abs(netOD).^(coef(3,:));
sigma_D=sqrt(netOD.^2.*sigma_coef(1,:).^2+(netOD.^2).^(coef(3,:)).*sigma_coef(2,:));

%% ELIMINACIÓN DE VARIABLES

clearvars -except coef alpha_cal beta_cal D_cal... 
                  Nx_I Ny_I...
                  OD0...
                  DD_pixel sigma_DD netODD NOD_D D_ref...
                  xmin xmax ymin ymax D_max...
                  Nx Ny N_data D_pixel sigma_D netOD

%% CÁLCULO DE LA MEDIA Y LA DESVIACIÓN DE LOS PARÁMETROS

% Magnitudes conocidas para la mancha de dosis conocida
dose=D_ref; % dosis del recorte de referencia
NOD=NOD_D;  % Calculada a partir de la dosis y la funcion de calibración

d_dose=coef(1,:)+coef(3,:).*coef(2,:).*NOD.^(coef(3,:)-1);  %derivada curva cal respecto NOD
d_mu_da=d_dose.*NOD;  %derivada de la serie de Taylor respecto alfa (a)
d_mu_db=d_dose.*OD0; %derivada de la serie de Taylor respecto beta (b)

% a = alfa b = beta i = coeficiente independiente
% DD_pixel dosis del pixel, sigma_DD sigma de dosis del pixel
Cai=sum((dose-DD_pixel).*d_mu_da./sigma_DD.^2,2); % suma sobre los 3 canales
Caa=sum(d_dose.*NOD.*d_mu_da./sigma_DD.^2,2);
Cab=sum(d_dose.*OD0.*d_mu_da./sigma_DD.^2,2);

Cbi=sum((dose-DD_pixel).*d_mu_db./sigma_DD.^2,2);
Cbb=sum(d_dose.*OD0.*d_mu_db./sigma_DD.^2,2);

% Resultado de alfa y beta del recorte de referencia
alpha=(Cai.*Cbb-Cbi.*Cab)./(Caa.*Cbb-Cab.^2);
beta=(Cai.*Cab-Cbi.*Caa)./(Cab.^2-Caa.*Cbb);

% Añadimos la información de alfa y beta a las alfas y betas calculadas en la calibración
N_limit=sum(D_cal<2)+1;
alpha=[reshape(alpha_cal(:,N_limit:end),1,[]), reshape(alpha,1,[])];
beta=[reshape(beta_cal(:,N_limit:end),1,[]), reshape(beta,1,[])];

% Calculamos las medias y sigmas
alpha_media=mean(alpha);
beta_media=mean(beta);
lambda_alpha=1./std(alpha).^2;
lambda_beta=1./std(beta).^2;

%% CÁLCULO DE LOS DIFERENTES PARÁMETROS (M. BROYDEN 3 VARIABLES)

% Inicialización de las variables
alpha=zeros(N_data,1);
beta=zeros(N_data,1);
NOD=zeros(N_data,3)+(netOD>0).*netOD+(netOD<0).*1e-5;

% Crear valor funcion y jacobiano
F=zeros(3,1); % Orden: NODR, NODG, NODB,
J=zeros(3); % Jacobiano

for mm=1:N_data  % N_data es el número de pixeles del area a la que aplicamos multicanal
    
    
    % CÁLCULO DEL JACOBIANO PARA INICIAR LA MATRIZ APROXIMADA DEL MÉTODO DE
    % BROYDEN.
    
    % Dosis a partir de NOD y sus derivadas respecto a la NOD
    dose=coef(1,:).*NOD(mm,:)+coef(2,:).*NOD(mm,:).^(coef(3,:));
    d_dose=coef(1,:)+coef(3,:).*coef(2,:).*NOD(mm,:).^(coef(3,:)-1); % derivada primera calibracion respecto NOD
    d2_dose=(coef(3,:)-1).*coef(3,:).*coef(2,:).*NOD(mm,:).^(coef(3,:)-2); % derivada segunda calibracion respecto NOD
    d3_dose=(coef(3,:)-2).*(coef(3,:)-1).*coef(3,:).*coef(2,:).*NOD(mm,:).^(coef(3,:)-3); % derivada tercera calibracion respecto NOD
        
    % Coeficientes cálculo de alpha y beta
    Ca=sum(d_dose.^2.*NOD(mm,:).^2./sigma_D(mm,:).^2)+lambda_alpha;
    Cb=sum(d_dose.^2.*OD0.^2./sigma_D(mm,:).^2)+lambda_beta;
    Cab=sum(d_dose.^2.*NOD(mm,:).*OD0./sigma_D(mm,:).^2);
    Cia=sum(d_dose.*NOD(mm,:).*(D_pixel(mm,:)-dose)./sigma_D(mm,:).^2)+lambda_alpha.*alpha_media;
    Cib=sum(d_dose.*OD0.*(D_pixel(mm,:)-dose)./sigma_D(mm,:).^2)+lambda_beta.*beta_media;
        
    % Derivadas de los coeficientes
    d_Ca=2.*d_dose.*NOD(mm,:).*(NOD(mm,:).*d2_dose+d_dose)./sigma_D(mm,:).^2;
    d_Cb=2.*d_dose.*d2_dose.*OD0.^2./sigma_D(mm,:).^2;
    d_Cab=d_dose.*OD0.*(2.*NOD(mm,:).*d2_dose+d_dose)./sigma_D(mm,:).^2;
    d_Cia=d2_dose.*NOD(mm,:).*(D_pixel(mm,:)-dose)./sigma_D(mm,:).^2+d_dose.*(D_pixel(mm,:)-dose)./sigma_D(mm,:).^2-d_dose.^2.*NOD(mm,:)./sigma_D(mm,:).^2;
    d_Cib=d2_dose.*OD0.*(D_pixel(mm,:)-dose)./sigma_D(mm,:).^2-d_dose.^2.*OD0./sigma_D(mm,:).^2;
        
    % Alpha y beta
    alpha(mm,1)=(Cia.*Cb-Cib.*Cab)./(Ca.*Cb-Cab.^2);
    beta(mm,1)=(Cia.*Cab-Cib.*Ca)./(Cab.^2-Ca.*Cb);
        
    % Derivadas de alpha y beta
    d_alpha=(d_Cia.*Cb+Cia.*d_Cb-d_Cib.*Cab-Cib.*d_Cab-alpha(mm,1).*(d_Ca.*Cb+Ca.*d_Cb-2.*Cab.*d_Cab))./(Ca.*Cb-Cab.^2);
    d_beta=(d_Cia.*Cab+Cia.*d_Cab-d_Cib.*Ca-Cib.*d_Ca-beta(mm,1).*(2.*Cab.*d_Cab-d_Ca.*Cb-Ca.*d_Cb))./(Cab.^2-Ca.*Cb);
        
    % Variación de la NOD respecto a la media
    var_NOD=NOD(mm,:).*alpha(mm,1)+OD0.*beta(mm,1);
    d_var_NOD=alpha(mm,1)+NOD(mm,:).*d_alpha+OD0.*d_beta;
        
    % Media estadística de la dosis  (Serie de Taylor)
    mu=dose+d_dose.*var_NOD;
    der_mu=d_dose.*(alpha(mm,1)+1)+d2_dose.*var_NOD;
    
    d_mu=d_dose.*(d_var_NOD+1)+d2_dose.*var_NOD; % derivada de mu respecto NOD
    d_der_mu=d2_dose.*(alpha(mm,1)+1+d_var_NOD)+d_dose.*d_alpha+d3_dose.*var_NOD; % derivada (para cálculo de Jacobiano) de la "derivada de mu respecto NOD (viene de minimización)" respecto NOD
        
    % Aplicamos método Newton (similar, Broyden)
	% Diferencia media dosis y dosis pesada por la incertidumbre al
    % cuadrado
    dif=(mu-D_pixel(mm,:))./sigma_D(mm,:).^2;
        
    % Cálculo del valor de la funcion (funciones correspondientes al sistema no lineal que nos queda por resolver)
    F(1)=sum(dif.*der_mu./d_dose);
    F(2)=dose(1,1)-dose(1,2);
    F(3)=dose(1,1)-dose(1,3);
                
    % Cálculo del jacobiano
    J(1,1:3)=d_mu.*der_mu./(d_dose.*sigma_D(mm,:).^2)+dif.*d_der_mu./d_dose-dif.*der_mu.*d2_dose./d_dose.^2; 
    J(2,1)=d_dose(1,1);
    J(2,2)=-d_dose(1,2);    
    J(3,1)=d_dose(1,1);
    J(3,3)=-d_dose(1,3);
    
    % Inversa de la matrix aproximada
    A=inv(J);
    
    % Inicialización bucle
    indice_max=0;
    dif_dosis=1;
    F_new=zeros(3,1);
    NOD_pix=NOD(mm,:);
    NOD_new=NOD_pix-(J\F)'; NOD_new(NOD_new<0)=1e-10;
    
    while indice_max<=50 && dif_dosis>1e-15
        
        %CÁLCULO DEL VALOR DE LA FUNCIÓN EN EL SIGUIENTE PASO
        
        % Dosis a partir de NOD y sus derivadas respecto a la NOD
        dose=coef(1,:).*NOD_new+coef(2,:).*NOD_new.^(coef(3,:));
        d_dose=coef(1,:)+coef(3,:).*coef(2,:).*NOD_new.^(coef(3,:)-1);
        d2_dose=(coef(3,:)-1).*coef(3,:).*coef(2,:).*NOD_new.^(coef(3,:)-2);
        
        % Coeficientes cálculo de alpha y beta
        Ca=sum(d_dose.^2.*NOD_new.^2./sigma_D(mm,:).^2)+lambda_alpha;
        Cb=sum(d_dose.^2.*OD0.^2./sigma_D(mm,:).^2)+lambda_beta;
        Cab=sum(d_dose.^2.*NOD_new.*OD0./sigma_D(mm,:).^2);
        Cia=sum(d_dose.*NOD_new.*(D_pixel(mm,:)-dose)./sigma_D(mm,:).^2)+lambda_alpha.*alpha_media;
        Cib=sum(d_dose.*OD0.*(D_pixel(mm,:)-dose)./sigma_D(mm,:).^2)+lambda_beta.*beta_media;
        
        % Alpha y beta
        alpha(mm,1)=(Cia.*Cb-Cib.*Cab)./(Ca.*Cb-Cab.^2);
        beta(mm,1)=(Cia.*Cab-Cib.*Ca)./(Cab.^2-Ca.*Cb);
        
        % Variación de la NOD respecto a la media
        var_NOD=NOD_new.*alpha(mm,1)+OD0.*beta(mm,1);
       
        % Media estadística de la dosis
        mu=dose+d_dose.*var_NOD;
        der_mu=d_dose.*(alpha(mm,1)+1)+d2_dose.*var_NOD;
        
        % Diferencia media dosis y dosis pesada por la incertidumbre al
        % cuadrado
        dif=(mu-D_pixel(mm,:))./sigma_D(mm,:).^2;
        
        % Cálculo del valor de la funcion
        F_new(1)=sum(dif.*der_mu./d_dose);
        F_new(2)=dose(1,1)-dose(1,2);
        F_new(3)=dose(1,1)-dose(1,3);
        
        
        %CÁLCULO DE LAS DIFERENCIAS
        dif_NOD=NOD_new-NOD_pix;
        dif_F=F_new-F;
        
        % Denominador actualización matriz aproximada
        denom=dif_NOD*A*dif_F;
        denom(denom==0)=1;
        
        %ACTUALIZACIÓN DE LA MATRIZ APROXIMADA
        A=A+(dif_NOD'-A*dif_F)*dif_NOD*A/denom;
        
        %ACTUALIZACIÓN DE LAS VARIABLES
        NOD_pix=NOD_new;
        NOD_new=NOD_pix-(A*F_new)'; NOD_new(NOD_new<0)=1e-10;
        
        F=F_new;
        dosis_new=coef(1,:).*NOD_new+coef(2,:).*NOD_new.^(coef(3,:));
        
        %ACTUALIZACIÓN CONDICIONES BUCLE
        dif_dosis=max(abs([dosis_new(1)-dosis_new(2),dosis_new(1)-dosis_new(3),dosis_new(2)-dosis_new(3)]))/min(dosis_new);
        indice_max=indice_max+1;
        
    end
   
    % GUARDAMOS EL VALOR DE LA NOD OBTENIDO
    NOD(mm,:)=NOD_new;
    
end

% CÁLCULO DE LA DOSIS
Dose=coef(1,:).*NOD+coef(2,:).*NOD.^(coef(3,:));
Dose(Dose>D_max)=D_max;

% Reordenación de los vectores en matrices
Dose_R=reshape(Dose(:,1),[Ny,Nx]);
Dose_G=reshape(Dose(:,2),[Ny,Nx]);
Dose_B=reshape(Dose(:,3),[Ny,Nx]);
D_R=reshape(D_pixel(:,1),[Ny,Nx]);
D_G=reshape(D_pixel(:,2),[Ny,Nx]);
D_B=reshape(D_pixel(:,3),Ny,Nx);

%% SAVE DATA

% Factor de escalamiento
f_sca=65535./D_max; % (niveles/Gy)

% Dosis sin multichannel
I_sin=zeros(Ny_I,Nx_I,3);

I_sin(ymin:ymax,xmin:xmax,1)=f_sca.*D_R;
I_sin(ymin:ymax,xmin:xmax,2)=f_sca.*D_G;
I_sin(ymin:ymax,xmin:xmax,3)=f_sca.*D_B;

I_sin=uint16(round(I_sin));

[archivo,ruta]=uiputfile('triple_canal.tif','Guardar mapa de dosis sin corregir');

% Pegamos el nombre del archivo y de la dirección
if archivo==0
 return;
else
 dat_archivo=strcat(ruta,archivo);
end

imwrite(I_sin,dat_archivo)

% Dosis multichannel
I_multi=zeros(Ny_I,Nx_I,3);

I_multi(ymin:ymax,xmin:xmax,1)=f_sca.*Dose_R;
I_multi(ymin:ymax,xmin:xmax,2)=f_sca.*Dose_G;
I_multi(ymin:ymax,xmin:xmax,3)=f_sca.*Dose_B;

I_multi=uint16(I_multi);

[archivo,ruta]=uiputfile('multicanal.tif','Guardar mapa de dosis corregido por el multicanal');

% Pegamos el nombre del archivo y de la dirección
if archivo==0
 return;
else
 dat_archivo=strcat(ruta,archivo);
end

imwrite(I_multi,dat_archivo)

end

