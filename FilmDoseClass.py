import ntpath
import numpy as np
from lmfit import Model
import tifffile #https://pypi.org/project/tifffile/
from PIL import Image
import FilmCalibrationClass_Dam
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import tkinter as tk
from tkinter import simpledialog
from tkinter import filedialog


class FilmDose:
    def __init__(self, filmfile):
        self.workingdir = ntpath.dirname(filmfile) + '\\'
        self.imagefilename = ntpath.basename(filmfile)
        image = tifffile.imread(filmfile)
        self.PILImage = Image.open(filmfile)
        tif = tifffile.TiffFile(filmfile)
        self.tif = tif
        # Get dpi and pixel size
        self.dpi = np.array([1.0, 1.0])
        tag = tif.pages[0].tags['XResolution']       
        self.dpi[1] = tag.value[0]/tag.value[1]
        tag = tif.pages[0].tags['YResolution']       
        self.dpi[0] = tag.value[0]/tag.value[1]
        if(str(tif.pages[0].tags['ResolutionUnit']).split()[-1]=='INCH'):
            self.pixelsize = 25.4/self.dpi
        elif(str(tif.pages[0].tags['ResolutionUnit']).split()[-1]=='MM'):
            self.dpi = self.dpi/25.4
            self.pixelsize = 1/self.dpi

        self.imagearray = np.array(image)
        self.OD = np.log10(65535.0/self.imagearray)
        self.OD0 = np.log10(np.array([65535.0, 65535.0, 65535.0]))
        self.ODref = np.log10(np.array([65535.0, 65535.0, 65535.0]))
        self.ROI = np.array(1)
        self.Dref = np.empty(3)
        self.CalibrationCorrectionFactors = np.array([1.0, 1.0, 1.0])

    def get_don(self):
        nod = np.clip(self.OD - self.OD0, 0.0)
        return nod

    def set_don(self, don):
        self.OD = don + self.OD0
        self.set_imagearray_from_od(self.OD)

    def get_do(self):
        return self.OD

    def set_do(self, do):
        self.OD = do
        self.set_imagearray_from_od(do)

    def get_dose(self):
        darray = np.zeros(self.OD.shape, dtype=np.float)
        #for c in [0, 1, 2]:
        #    darray[:, :, c] = self.Calibration.get_dose_from_nod((self.OD[:, :, c] - self.OD0[c]), c)
        darray = self.Calibration.get_dose_from_nod((self.OD - self.OD0))
        # Se permiten dosis hasta un 20% superiores al máximo de la curva de calibración
        darray = np.clip(darray, 0, 1.2 * np.max(self.Calibration.D))
        return darray

    def get_pixel(self):
        return self.imagearray

    def set_pixel(self, pixel):
        self.imagearray = pixel
        self.set_od_from_imagearray(pixel)


    def set_od_from_imagearray(self, imarr):
        self.imagearray = np.array(imarr)
        self.OD = np.log10(65535.0 / self.imagearray)
        self.ROI = np.array(1)

    def set_imagearray_from_od(self, od):
        self.imagearray = np.clip(65535.0 / np.power(10, self.OD), 0, 65535).astype(np.uint16)

    def get_pixelheight(self):
        return self.pixelsize[0]
    
    def get_pixelwidth(self):
        return self.pixelsize[1]
       
    def set_pixelsize(self, ps):
        self.pixelsize[0] = ps
        self.pixelsize[1] = ps
        self.dpi = 25.4 / self.pixelsize
        
    def set_pixelheight(self,ph):
        self.pixelsize[0] = ph
        self.dpi[0] = 25.4 / self.pixelsize[0]
    def set_pixelwidth(self,pw):
        self.pixelsize[1] = pw   
        self.dpi[1] = 25.4 / self.pixelsize[1]

    # Genera la matriz de dosis a partir de la matriz de densidad óptica y la calibración
    def DoseArrays(self):
        darray = np.zeros(self.OD.shape, dtype=np.float)
        for c in [0, 1, 2]:
            darray[:, :, c] = self.Calibration.get_dose_from_nod((self.OD[:, :, c] - self.OD0[c]), c)
        # Se permiten dosis hasta un 20% superiores al máximo de la curva de calibración
        darray = np.clip(darray, 0, 1.2 * np.max(self.Calibration.D))
        return darray

    # Guarda la imagen en formato tiff 48-bits con valor de pixel lineal con la dosis
    def SaveDoseImage(self, outputfile, Dmax):
        arr = (self.DoseArrays() * 65535 / Dmax)
        arr[arr > 65535] = 65535
        darray = arr.astype(np.uint16)
        tifffile.imwrite(outputfile, darray, resolution=(self.dpi[0], self.dpi[1]))

    # Guarda la imagen corregida en Tiff 48-bits
    def SaveCorrectedPVImage(self, outputfile):
        PV = (65535.0 / np.power(10, self.OD)).astype(np.uint16)
        tifffile.imwrite(outputfile, PV, resolution=(self.dpi[0], self.dpi[1]))

    def show_Image(self):
        fig = plt.figure()
        plt_image = plt.imshow(self.PILImage)
        plt.show()
    
    def show_DoseImage(self, dmax):
        fig = plt.figure()
        plt_image = plt.imshow(self.generate_DoseImage(dmax))
        plt.show()

    def generate_DoseImage(self, dmax):
        arr = np.clip(self.DoseArrays()*256/dmax,0,255).astype(np.uint8)
        return arr
    
    def add_calibration_objects(self, fcalobj):
        self.Calibration = fcalobj
        self.OD0 = self.Calibration.OD0

    def add_calibration(self, pv, d, pv_std_dev=None):
        fcalobj = FilmCalibrationClass_Dam.FilmCalibration()
        fcalobj.calibration_from_arrays(pv, d, pv_std_dev)
        self.Calibration = fcalobj
        self.OD0 = self.Calibration.OD0

    # Carga una curva de calibración desde un fichero
    def AddCalibrationFromFile(self, inputfile):
        #root = tk.Tk()
        #root.withdraw()
        #file_path = filedialog.askopenfilename(title='Open Calibration file')
        #root.destroy()
        self.Calibration = FilmCalibrationClass_Dam.FilmCalibration(inputfile)
        self.OD0 = self.Calibration.OD0
        
    # Función que permite seleccionar una ROI rectangular en la imagen.
    # Devuelve las coordenadas de la ROI
    def SelectRectangle(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt_image = plt.imshow(self.PILImage)
        rs = widgets.RectangleSelector(
            ax, self.__onselect, drawtype='box',
            rectprops=dict(facecolor='red', edgecolor='black', alpha=0.2, fill=True))
        plt.show(block=True)
        x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
        y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
        rect = np.array([x, y])
        return rs
        
    # Función que captura las coordenadas del mouse al hacer click en una figura.
    # Se pasa como argumento a SelectRectangle
    def __onselect(self,eclick, erelease):
        if eclick.ydata > erelease.ydata:
            eclick.ydata, erelease.ydata=erelease.ydata, eclick.ydata
        if eclick.xdata > erelease.xdata:
            eclick.xdata, erelease.xdata=erelease.xdata, eclick.xdata
        plt.close()
        return 
        
    # Solicita seleccionar ROIs de un área no irradiada y de otra irradiada con dosis conocida.
    # Aplica un factor de corrección a la calibración
    def CalibrationCorrection(self):
        import ctypes  # An included library with Python install.   

        # Pide seleccionar el recorte sin irradiar
        ctypes.windll.user32.MessageBoxW(0, "Select non irradiated ROI", "", 0)
        rs = self.SelectRectangle()
        x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
        y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
        self.ROI = self.OD[y[0]:(y[1]+1), x[0]:(x[1]+1), 0:3]
        od0 = np.array(self.ROI)

        # Pide seleccionar el recorte irradiado a dosis conocida de referencia
        ctypes.windll.user32.MessageBoxW(0, "Select irradiated ROI", "", 0)
        rs = self.SelectRectangle()
        x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
        y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
        #self.ROI = self.imagearray[y[0]:(y[1]+1),x[0]:(x[1]+1),0:3]
        #self.ODref = np.mean(np.log10(65535/self.ROI),axis=(0,1))
        self.ROI = self.OD[y[0]:(y[1]+1),x[0]:(x[1]+1),0:3]
        odref = np.array(self.ROI)

        # Solicita el valor de la dosis de referencia
        dref = np.empty(3)
        application_window = tk.Tk()
        application_window.withdraw()
        answer = simpledialog.askstring("Input", "Introduce de reference dose (Gy)",
                                        parent=application_window)
        dref.fill(np.float(answer))
        application_window.destroy()
        # Calcula la ODnet que corresponde a la dosis de referencia usando la curva de calibración.
        self.calculate_calibration_correction_factors(od0, odref, dref)
        return True

    def calculate_calibration_correction_factors(self, od0, odref, dref):
        # Calcula la ODnet que corresponde a la dosis de referencia usando la curva de calibración.
        self.Dref = dref
        self.OD0 = np.mean(od0,axis=(0,1))
        self.ODref = np.mean(odref,axis=(0,1))
        ODnetRefCal = np.array([0.0, 0.0, 0.0])
        for c in np.array([0, 1, 2]):
            ODnetRefCal[c] = self.Calibration.GetODnetFromDose(self.Dref[c], c)

        # Calcula el factor de corrección de dosis
        self.CalibrationCorrectionFactors = ODnetRefCal / (self.ODref - self.OD0)
        # Aplica el factor de corrección. Llevamos la película a las condiciones de la calibración
        self.OD = self.Calibration.OD0 + (self.OD - self.OD0) * self.CalibrationCorrectionFactors
        self.OD0 = self.Calibration.OD0

        # Calcula las alphas y betas del roi de referencia para dosimetria multicanal
        coef1 = self.Calibration.DevicParam_A
        coef2 = self.Calibration.DevicParam_B
        coef3 = self.Calibration.DevicParam_n
        sigma_coef1 = self.Calibration.Sigma_A
        sigma_coef2 = self.Calibration.Sigma_B

        D = dref[0]
        nod_ref = odref-np.mean(od0,axis=(0,1))
        NOD = np.array([self.Calibration.GetODnetFromDose(D, c) for c in [0, 1, 2]])
        Dcal_pixel = coef1 * nod_ref + \
                     np.sign(nod_ref) * coef2 * np.power(np.abs(nod_ref), coef3)
        sigma_Dcal = np.sqrt(np.power(np.abs(nod_ref) * sigma_coef1, 2) +
                             np.power(np.abs(nod_ref), 2 * coef3) * np.power(sigma_coef2, 2))

        d_dose = coef1 + coef3 * coef2 * np.power(NOD, (coef3 - 1))
        d_mu_da = d_dose * NOD
        d_mu_db = d_dose * self.OD0

        Cai = np.sum((D - Dcal_pixel) * d_mu_da / np.power(sigma_Dcal, 2), axis=2)
        Caa = np.sum(d_dose * NOD * d_mu_da / np.power(sigma_Dcal, 2), axis=2)
        Cab = np.sum(d_dose * self.OD0 * d_mu_da / np.power(sigma_Dcal, 2), axis=2)

        Cbi = np.sum((D - Dcal_pixel) * d_mu_db / np.power(sigma_Dcal, 2), axis=2)
        Cbb = np.sum(d_dose * self.OD0 * d_mu_db / np.power(sigma_Dcal, 2), axis=2)

        alpha = np.mean((Cai * Cbb - Cbi * Cab) / (Caa * Cbb - Cab * Cab))
        sigma_alpha = np.std((Cai * Cbb - Cbi * Cab) / (Caa * Cbb - Cab * Cab))
        beta = np.mean((Cai * Cab - Cbi * Caa) / (Cab * Cab - Caa * Cbb))
        sigma_beta = np.std((Cai * Cab - Cbi * Caa) / (Cab * Cab - Caa * Cbb))
        npoints = (nod_ref.shape[0] * nod_ref.shape[1])
        self.Calibration.AlphaCal = (self.Calibration.AlphaCal*self.Calibration.AlphaBetaPoints +
                                     alpha*npoints)/(npoints + self.Calibration.AlphaBetaPoints)
        self.Calibration.BetaCal = (self.Calibration.BetaCal * self.Calibration.AlphaBetaPoints +
                                     beta * npoints) / (npoints + self.Calibration.AlphaBetaPoints)

        self.Calibration.SigmaAlphaCal = np.sqrt((np.power(self.Calibration.SigmaAlphaCal, 2) *
                                                  self.Calibration.AlphaBetaPoints +
                                                  sigma_alpha*sigma_alpha * npoints) /
                                                 (npoints + self.Calibration.AlphaBetaPoints))
        self.Calibration.SigmaBetaCal = np.sqrt((np.power(self.Calibration.SigmaBetaCal, 2) *
                                                  self.Calibration.AlphaBetaPoints +
                                                  sigma_beta * sigma_beta * npoints) /
                                                 (npoints + self.Calibration.AlphaBetaPoints))
        if self.Calibration.AlphaBetaPoints == 0:
            self.Calibration.AlphaCal = 0.0
            self.Calibration.BetaCal = 0.0

    # Carga los ficheros con los polinomios de corrección de flatscanner y genera una lista con los coeficientes
    # Se carga un fichero por cada nivel de dosis, en orden creciente con la dosis.
    def LoadFlatScanCorrFiles(self):
        self.FlatScanCoeffList = []
        count = 0
        file_path = 'something'
        while True:
            appwindow = tk.Tk()
            appwindow.withdraw()
            file_path = filedialog.askopenfilename(title='Select FlatScan file {}.'.format(count + 1))
            appwindow.destroy()
            count = count + 1
            if not file_path:
                break
            coeff = np.array(0.001*np.loadtxt(file_path, usecols=(0, 1, 2), skiprows=0))
            self.FlatScanCoeffList.append(coeff)

    # Aplica la corrección de flatscanner.
    # Sólo se aplica en la dirección transversal al movimiento de la lámpara
    def ApplyFlatScanCorr(self):
        self.FScorr = np.ones(self.imagearray.shape)
        nlevels = len(self.FlatScanCoeffList)  # numero de niveles de dosis
        odcenter = np.ones([len(self.FlatScanCoeffList), 3])
        center = self.imagearray.shape[0]*self.pixelsize[0]*0.5  # coordenada del centro de la imagen en altura
        # calcula densidad óptica en el centro
        for i in range(len(self.FlatScanCoeffList)):
            odcenter[i,:] = np.polyval(self.FlatScanCoeffList[i], center)
        coeffarray=np.ones([len(self.FlatScanCoeffList),
                            self.FlatScanCoeffList[0].shape[0], 3])
        # Pasa la lista de coeficientes a un array
        for i in np.arange(coeffarray.shape[0]):
            for j in np.arange(coeffarray.shape[1]):
                for c in np.arange(3):
                    coeffarray[i, j, c] = self.FlatScanCoeffList[i][j, c]
        # Inicializa el array que va a contener los factores de corrección para cada píxel
        farray = np.ones([nlevels, self.imagearray.shape[0], 3])
        for h in np.arange(self.imagearray.shape[0]):
            ypos = h*self.pixelsize[0]
            for c in np.arange(3):
                for nl in np.arange(nlevels):
                    # Factores de corrección para cada nivel de dosis
                    farray[nl,h,c] = np.polyval(coeffarray[nl,:,c],ypos)             
                for w in np.arange(self.imagearray.shape[1]):
                    # En cada pixel, interpola para obtener el factor de corrección correspondiente a su
                    # densidad optica
                    fcorr = self.__computeFScorrfactor(self.OD[h, w, c],
                                                        ypos,
                                                        farray[:, h, c],
                                                        odcenter[:, c])
                    self.FScorr[h, w, c] = fcorr
                    self.OD[h, w, c] = self.OD[h, w, c]/fcorr  # aplica la corrección

        # Genera una imagen tif 24-bits con el mapa de corrección obtenido
        FSimage = np.array((255*(self.FScorr-0.9)/0.2).astype(np.uint8))
        imFS = Image.fromarray(FSimage)
        imname = 'ScanCorrMap_' + self.imagefilename
        imFS.save(self.workingdir+imname)           
        return self.FScorr
    
    def __computeFScorrfactor(self, od, x, farr, odCentral):
        index = 0
        fcorr = 1.0
        nlevels = farr.shape[0]
        for i in np.arange(nlevels):
            odUnif = farr[i]
            if od > odUnif:
                index = i
        if index == 0 or index == (nlevels-1):
            fcorr = farr[index]/odCentral[index]
        if 0 < index < (nlevels - 1):
            od1 = farr[index]
            od2 = farr[index+1]
            f1 = od1/odCentral[index]
            f2 = od2/odCentral[index+1]
            a = (f2-f1)/(od2-od1)
            b = f1 - a*od1
            fcorr = a*od + b
        return fcorr

    # Método de corrección multicanal para un parámetro multiplicativo a la densidad óptica.
    # Equivalente al implementado en el programa de matlab.
    def multichannel_correction1(self):
        from scipy import optimize
        import ctypes  # An included library with Python install.
        # Selecciona el área a la que se va a aplicar la corrección
        ctypes.windll.user32.MessageBoxW(0, "Select area for multichannel correction", "", 0)
        rs = self.SelectRectangle()
        x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
        y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
        od_roi = self.OD[y[0]:(y[1] + 1), x[0]:(x[1] + 1), 0:3]
        self.alpha = np.ones([self.OD.shape[0], self.OD.shape[1]])
        lim = 0.2
        cont = 0.0
        cont2 = 0.0
        cont_lim = od_roi.shape[0]*od_roi.shape[1]*od_roi.shape[2]
        for h in np.arange(od_roi.shape[0]):
            if cont2/cont_lim > 0.05:
                cont2 = 0.0
                print(f'Multichannel correction process: {np.trunc(100*cont/cont_lim)}%')
            for w in np.arange(od_roi.shape[1]):
                #fittet_params = optimize.minimize(self.optimization_func1,
                #                                  x0=np.array([1.0]),
                #                                  args=(od_roi[h, w, :], self.OD0[:])).x
                fittet_params = self.mi_optimizacion(np.array([1.0]), od_roi[h, w, :],self.OD0[:])
                if fittet_params[0] > 1.2:
                    fittet_params[0] = 1.2
                if fittet_params[0] < 0.8:
                    fittet_params[0] = 0.8
                self.alpha[h + y[0], w + x[0]] = fittet_params[0]
                for c in np.arange(3):
                    od = od_roi[h, w, c]*fittet_params[0]
                    self.OD[h + y[0], w + x[0], c] = od
                    cont = cont + 1.0
                    cont2 = cont2 + 1.0
        # Genera una imagen con el mapa de corrección multicanal obtenido
        alpha_image = np.array((65535 * 0.5 + self.alpha * 1000))
        np.clip(alpha_image, 0, 65535)
        imname = 'MultichanelMap_' + self.imagefilename
        tifffile.imwrite(self.workingdir + imname,
                         alpha_image.astype(np.uint16),
                         resolution=(self.dpi[0], self.dpi[1]))

    def mi_optimizacion(self, params, od, od0):
        paso = 0.002
        #Determina dirección
        if self.optimization_func1(params+paso, od, od0) > self.optimization_func1(params-paso, od, od0):
            paso = -1*paso
        nparams = params
        f1 = self.optimization_func1(nparams, od, od0)
        f2 = self.optimization_func1(nparams+paso, od, od0)
        while f2 < f1:
            f1 = f2
            nparams = nparams + paso
            if abs(nparams[0] - 1.0) > 0.1:
                break
            f2 = self.optimization_func1(nparams + paso, od, od0)
        return nparams

    # Función que se optimiza en el método de corrección multicanal1
    def optimization_func1(self, params, od, od0):
        nod = od*params[0]-od0
        nod[nod < 0] = 0.0
        if (1-params[0]) > 0.2:
            return 10000.0
        dr, dg, db = self.Calibration.get_dose_from_nod(nod)
        f = 3*np.power(dr - dg, 2) + np.power(dr - db, 2) + np.power(dg - db, 2)
        return f

    # Método de corrección multicanal con dos parámetros (alfa, beta).
    # Basado en el método desarrollado por Damian.
    # En desarrollo...
    def multichannel_correction_Dam(self):
        from scipy import optimize
        import ctypes  # An included library with Python install.
        ctypes.windll.user32.MessageBoxW(0, "Select area for multichannel correction", "", 0)
        rs = self.SelectRectangle()
        x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
        y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
        #x[0] = 650
        #x[1]=7104
        #y[0] = 225
        #y[1] = 293
        nod_roi = self.OD[y[0]:(y[1] + 1), x[0]:(x[1] + 1), 0:3]
        nod_roi = nod_roi - self.OD0
        self.alpha = np.zeros([self.OD.shape[0], self.OD.shape[1]])
        self.beta = np.zeros([self.OD.shape[0], self.OD.shape[1]])
        lim = 0.2
        cont = 0.0
        cont2 = 0.0
        cont_lim = nod_roi.shape[0] * nod_roi.shape[1] * nod_roi.shape[2]
        F=np.zeros([3, 1])
        J = np.zeros([3,3])

        coef1 = self.Calibration.DevicParam_A
        coef2 = self.Calibration.DevicParam_B
        coef3 = self.Calibration.DevicParam_n
        alpha_media = self.Calibration.AlphaCal
        beta_media = self.Calibration.BetaCal
        lambda_alpha = 1.0 / np.power(self.Calibration.SigmaAlphaCal,2)
        lambda_beta = 1.0 / np.power(self.Calibration.SigmaBetaCal,2)
        sigma_coef1 = self.Calibration.Sigma_A
        sigma_coef2 = self.Calibration.Sigma_B
        sigma_coef3 = self.Calibration.Sigma_n

        od0 = self.OD0

        indice_medio = 0
        alpha_average = 0
        beta_average = 0

        for h in range(nod_roi.shape[0]):
            if cont2/cont_lim > 0.05:
                cont2 = 0.0
                print(f'Multichannel correction process: {np.trunc(100*cont/cont_lim)}%')
            for w in range(nod_roi.shape[1]):
                netOD = np.array([nod_roi[h, w, c] for c in [0, 1, 2]])
                D_pixel = coef1*netOD+np.sign(netOD)*coef2*np.power(np.abs(netOD),coef3)
                sigma_D = np.sqrt(np.power(netOD*sigma_coef1, 2) +
                                  np.power(np.abs(netOD), 2*coef3)*np.power(sigma_coef2, 2))
                nod = np.array([0.0 + (netOD[c] > 0) * netOD[c] + (netOD[c] < 0) * 0.00001 for c in [0, 1, 2]])
                dose = coef1*nod+coef2*np.power(nod,coef3)
                d_dose = coef1 + coef3*coef2*np.power(nod,coef3-1)
                d2_dose = (coef3-1)*coef3*coef2*np.power(nod,coef3-2)
                d3_dose = (coef3 - 2)*(coef3 - 1) * coef3 * coef2 * np.power(nod, coef3 - 3)

                Ca = np.sum(np.power(d_dose * nod / sigma_D, 2)) + lambda_alpha
                Cb = np.sum(np.power(d_dose * od0 / sigma_D, 2)) + lambda_beta
                Cab = np.sum(np.power(d_dose / sigma_D, 2) * nod * od0)
                Cia = np.sum(d_dose * nod * (D_pixel - dose) / np.power(sigma_D, 2)) + lambda_alpha * alpha_media
                Cib = np.sum(d_dose * od0 * (D_pixel - dose) / np.power(sigma_D, 2)) + lambda_beta * beta_media

                d_Ca = 2 * d_dose * nod * (nod * d2_dose+d_dose) / np.power(sigma_D,2)
                d_Cb = 2 * d_dose * d2_dose * np.power(od0/sigma_D, 2)
                d_Cab = d_dose * od0 * (2 * nod * d2_dose + d_dose) / np.power(sigma_D, 2)
                d_Cia = (d2_dose * nod*(D_pixel - dose) + d_dose * (D_pixel - dose) -
                        np.power(d_dose, 2) * nod) / np.power(sigma_D, 2)
                d_Cib = (d2_dose * od0 * (D_pixel - dose) - np.power(d_dose, 2) * od0) / np.power(sigma_D, 2)

                alpha = (Cia * Cb - Cib * Cab) / (Ca * Cb - Cab * Cab)
                beta = (Cia * Cab - Cib * Ca) / (Cab * Cab - Ca * Cb)

                if np.isnan(alpha) or np.isnan(beta):
                    xx=1

                d_alpha = (d_Cia * Cb + Cia * d_Cb - d_Cib * Cab - Cib * d_Cab -
                           alpha * (d_Ca * Cb + Ca * d_Cb - 2. * Cab * d_Cab)) / (Ca * Cb - Cab * Cab)
                d_beta = (d_Cia * Cab + Cia * d_Cab - d_Cib * Ca - Cib * d_Ca -
                          beta * (2. * Cab * d_Cab - d_Ca * Cb - Ca * d_Cb)) / (Cab * Cab - Ca * Cb)

                var_NOD = nod*alpha + od0 * beta
                d_var_NOD = alpha + nod * d_alpha + od0 * d_beta

                mu = dose + d_dose * var_NOD
                der_mu = d_dose * (alpha + 1) + d2_dose * var_NOD

                d_mu = d_dose * (d_var_NOD + 1) + d2_dose * var_NOD
                d_der_mu = d2_dose * (alpha + 1 + d_var_NOD) + d_dose * d_alpha + d3_dose * var_NOD

                dif = (mu - D_pixel) / np.power(sigma_D, 2)

                F[0, 0] = np.sum(dif * der_mu / d_dose)
                F[1, 0] = dose[0] - dose[1]
                F[2, 0] = dose[0] - dose[2]

                for c in np.arange(3):
                    J[0, c] = d_mu[c] * der_mu[c] / (d_dose[c] * np.power(sigma_D[c], 2))+dif[c] * d_der_mu[c] / \
                            d_dose[c] - dif[c] * der_mu[c] * d2_dose[c] / np.power(d_dose[c], 2)
                J[1, 0] = d_dose[0]
                J[1, 1] = -d_dose[1]
                J[2, 0] = d_dose[0]
                J[2, 2] = -d_dose[2]

                A = np.linalg.inv(J)

                indice_max = 0
                dif_dosis = 1
                F_new = np.zeros([3, 1])
                NOD_pix = np.array([0.0, 0.0, 0.0]) +  nod
                NOD_new = NOD_pix - np.matrix.transpose(np.linalg.solve(J, F))
                NOD_new[NOD_new < 0] = 1e-10

                while (indice_max <= 30) and (dif_dosis > 1.0e-10):
                    nod = np.array([0.0, 0.0, 0.0]) + NOD_new
                    dose = coef1 * nod + coef2 * np.power(nod, coef3)
                    d_dose = coef1 + coef3 * coef2 * np.power(nod, coef3 - 1)
                    d2_dose = (coef3 - 1) * coef3 * coef2 * np.power(nod, coef3 - 2)

                    Ca = np.sum(np.power(d_dose * nod / sigma_D, 2), axis=1)[0] + lambda_alpha
                    Cb = np.sum(np.power(d_dose * od0 / sigma_D, 2),axis=1)[0] + lambda_beta
                    Cab = np.sum(np.power(d_dose / sigma_D, 2) * nod * od0,axis=1)[0]
                    Cia = np.sum(d_dose * nod * (D_pixel - dose) / np.power(sigma_D, 2),axis=1)[0] + \
                          lambda_alpha * alpha_media
                    Cib = np.sum(d_dose * od0 * (D_pixel - dose) / np.power(sigma_D, 2),axis=1)[0] + \
                          lambda_beta * beta_media

                    alpha = (Cia * Cb - Cib * Cab) / (Ca * Cb - Cab * Cab)
                    beta = (Cia * Cab - Cib * Ca) / (Cab * Cab - Ca * Cb)

                    var_NOD = nod * alpha + od0 * beta

                    mu = dose + d_dose * var_NOD
                    der_mu = d_dose * (alpha + 1) + d2_dose * var_NOD

                    dif = (mu - D_pixel) / np.power(sigma_D, 2)

                    F_new[0,0] = np.sum(dif * der_mu / d_dose, axis=1)[0]
                    F_new[1,0] = dose[0,0] - dose[0,1]
                    F_new[2,0] = dose[0,0] - dose[0,2]

                    dif_NOD = NOD_new - NOD_pix
                    dif_F = F_new - F

                    denom = np.matmul(dif_NOD, np.matmul(A,dif_F))[0]
                    if denom == 0: denom = 1

                    #dif_NODt = np.zeros([1,3])
                    #dif_NODt[0,:] = dif_NOD[:]
                    A = A + np.matmul(np.matmul(np.transpose(dif_NOD)-np.matmul(A,dif_F),dif_NOD), A)/denom

                    NOD_pix = np.array([0.0, 0.0, 0.0]) + NOD_new
                    NOD_new = NOD_pix - np.matrix.transpose(np.matmul(A, F_new))
                    NOD_new[NOD_new < 0] = 1e-10

                    F = np.zeros([3, 1]) + F_new
                    dosis_new = coef1 * NOD_new + coef2 * np.power(NOD_new, coef3)

                    dif_dosis = np.max(np.abs(np.array([dosis_new[0,0] - dosis_new[0,1],
                                                        dosis_new[0,0] - dosis_new[0,2],
                                                        dosis_new[0,1] - dosis_new[0,2]]))) / np.min(dosis_new)
                    indice_max = indice_max + 1

                indice_medio = indice_medio + indice_max
                alpha_average = alpha_average + alpha
                beta_average = beta_average + beta
                self.alpha[h + y[0], w + x[0]] = alpha
                self.beta[h + y[0], w + x[0]] = beta

                for c in [0,1,2]:
                    #aux = NOD_new[0, c] + od0[c]
                    self.OD[h + y[0], w + x[0], c] = NOD_new[0, c] + od0[c]
                    #self.OD[h + y[0], w + x[0], c] = (self.OD[h + y[0], w + x[0], c]-
                    #                                 od0[c]*(1+beta))/(1.0+alpha) + od0[c]
                cont2 = cont2 + 3.0
                cont = cont + 3.0

        indice_medio = indice_medio/(nod_roi.shape[0]*nod_roi.shape[1])
        alpha_average = alpha_average/(nod_roi.shape[0]*nod_roi.shape[1])
        beta_average = beta_average/(nod_roi.shape[0]*nod_roi.shape[1])
        print(f'alpha average = {alpha_average}')
        print(f'beta average = {beta_average}')
        print(f'indice medio = {indice_medio}')
        alpha_image = np.array((65535 * 0.5 + self.alpha*1000))
        np.clip(alpha_image, 0, 65535)
        imname = 'AlphaMap_' + self.imagefilename
        tifffile.imwrite(self.workingdir + imname,
                         alpha_image.astype(np.uint16),
                         resolution=(self.dpi[0], self.dpi[1]))

        beta_image = np.array((65535 * 0.5 + self.beta * 1000))
        np.clip(beta_image, 0, 65535)
        imname = 'BetaMap_' + self.imagefilename
        tifffile.imwrite(self.workingdir + imname,
                         beta_image.astype(np.uint16),
                         resolution=(self.dpi[0], self.dpi[1]))


    def optimization_func2a(self, params, odnet):

        nod = (odnet - params[1] * self.OD0)/(1 + params[0])
        dr = self.Calibration.get_dose_from_nod(nod[0], 0)
        dg = self.Calibration.get_dose_from_nod(nod[1], 1)
        db = self.Calibration.get_dose_from_nod(nod[2], 2)
        f = np.power(dr - dg, 2) + np.power(dr - db, 2) + np.power(params[0], 2) + np.power(params[1], 2)
        return f

    def optimization_func2b(self, params, odnet):

        nod = (odnet - params[1] * self.OD0)/(1 + params[0])
        dk = np.array([self.Calibration.get_dose_from_nod(odnet[c], c) for c in [0, 1, 2]])
        d = np.array([self.Calibration.get_dose_from_nod(nod[c], c) for c in [0, 1, 2]])
        delta_d = np.array([self.Calibration.Get_DerivateFromODnet(nod[c], c) for c in [0, 1, 2]])
        delta_d = delta_d*(params[0]*odnet + params[1]*self.OD0)/(1 + params[0])
        f = 0.0
        unc = np.array([1.0, 1.0, 10.0])
        for c in [0, 1, 2]:
            f = f + np.power(d[c] + delta_d[c] - dk[c], 2)/unc[c]
        #f = f + 0.001*np.power(params[0], 2) + 0.001*np.power(params[1], 2)
        f = f + np.power(d[0]-d[1], 2) + 0.1*np.power(d[0]-d[2], 2)
        return f
        

    