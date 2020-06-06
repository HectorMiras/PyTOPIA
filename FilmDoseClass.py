import ntpath
import numpy as np
from lmfit import Model
import tifffile #https://pypi.org/project/tifffile/
from PIL import Image
import FilmCalibrationClass
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
        

    
    def get_pixelheight(self):
        return self.pixelsize[0]
    
    def get_pixelwidth(self):
        return self.pixelsize[1]
       
    def set_pixelsize(self,ps):
        self.pixelsize[0] = ps
        self.pixelsize[1] = ps
        self.dpi = 25.4 / self.pixelsize
        
    def set_pixelheight(self,ph):
        self.pixelsize[0] = ph
        self.dpi[0] = 25.4 / self.pixelsize[0]
    def set_pixelwidth(self,pw):
        self.pixelsize[1] = pw   
        self.dpi[1] = 25.4 / self.pixelsize[1]
        
    def get_ImageRGB_24bits(self):
        return self.PILImage

    def get_nod(self):
        nod = self.OD-self.OD0
        nod[nod < 0] = 0.0
        return nod
    
    def show_Image(self):
        fig = plt.figure()
        plt_image = plt.imshow(self.PILImage)
        plt.show()
    
    def show_DoseImage(self, dmax):
        fig = plt.figure()
        plt_image = plt.imshow(self.generate_DoseImage(dmax))
        plt.show()

    def generate_DoseImage(self, dmax):
        arr = (self.DoseArrays()*256/dmax)
        arr[arr > 255] = 255
        arr = arr.astype(np.uint8)
        return arr
    
    # Carga una curva de calibración desde un fichero
    def AddCalibrationFromFile(self, inputfile):
        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename(title='Open Calibration file')
        root.destroy()
        self.Calibration = FilmCalibrationClass.FilmCalibration(file_path)
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
        #plt.pause(5)
        plt.show(block=True)  # Me bloquea el programa aunque cierre la ventana
        #plt.show()
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
        self.OD0 = np.mean(self.ROI, axis=(0, 1))

        # Pide seleccionar el recorte irradiado a dosis conocida de referencia
        ctypes.windll.user32.MessageBoxW(0, "Select irradiated ROI", "", 0)
        rs = self.SelectRectangle()
        x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
        y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
        #self.ROI = self.imagearray[y[0]:(y[1]+1),x[0]:(x[1]+1),0:3]
        #self.ODref = np.mean(np.log10(65535/self.ROI),axis=(0,1))
        self.ROI = self.OD[y[0]:(y[1]+1),x[0]:(x[1]+1),0:3]
        self.ODref = np.mean(self.ROI,axis=(0,1))

        # Solicita el valor de la dosis de referencia
        application_window = tk.Tk()
        application_window.withdraw()
        answer = simpledialog.askstring("Input", "Introduce de reference dose (Gy)",
                                        parent=application_window)
        self.Dref.fill(np.float(answer))
        application_window.destroy()
        # Calcula la ODnet que corresponde a la dosis de referencia usando la curva de calibración.
        ODnetRefCal = np.array([0.0,0.0,0.0])
        for c in np.array([0,1,2]): 
            ODnetRefCal[c] = self.Calibration.GetODnetFromDose(self.Dref[c],c)

        # Calcula el factor de corrección de dosis
        self.CalibrationCorrectionFactors = ODnetRefCal/(self.ODref-self.OD0)
        # Aplica el factor de corrección. Llevamos la película a las condiciones de la calibración
        self.OD = self.Calibration.OD0 + (self.OD-self.OD0)*self.CalibrationCorrectionFactors
        self.OD0 = self.Calibration.OD0
        return True

    # Genera la matriz de dosis a partir de la matriz de densidad óptica y la calibración
    def DoseArrays(self):
        darray = np.zeros(self.OD.shape, dtype=np.float)
        for c in [0, 1, 2]:
            darray[:, :, c] = self.Calibration.Get_DoseFromODnet((self.OD[:, :, c]-self.OD0[c]), c)
        # Se permiten dosis hasta un 20% superiores al máximo de la curva de calibración
        darray[darray > 1.2*np.max(self.Calibration.D)] = 1.2*np.max(self.Calibration.D)
        return darray
    
    # Guarda la imagen en formato tiff 48-bits con valor de pixel lineal con la dosis
    def SaveDoseImage(self, outputfile, Dmax):
        arr = (self.DoseArrays()*65535/Dmax)
        arr[arr > 65535] = 65535
        darray = arr.astype(np.uint16)
        tifffile.imwrite(outputfile, darray, resolution=(self.dpi[0], self.dpi[1]))

    # Guarda la imagen corregida en Tiff 48-bits
    def SaveCorrectedPVImage(self,outputfile):
        PV = (65535.0/np.power(10,self.OD)).astype(np.uint16)
        tifffile.imwrite(outputfile, PV, resolution=(self.dpi[0], self.dpi[1]))
        
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
        for i in range(coeffarray.shape[0]):
            for j in range(coeffarray.shape[1]):
                for c in range(3):
                    coeffarray[i, j, c] = self.FlatScanCoeffList[i][j, c]
        # Inicializa el array que va a contener los factores de corrección para cada píxel
        farray = np.ones([nlevels, self.imagearray.shape[0], 3])
        for h in range(self.imagearray.shape[0]):
            ypos = h*self.pixelsize[0]
            for c in [0,1,2]: 
                for nl in range(nlevels):
                    # Factores de corrección para cada nivel de dosis
                    farray[nl,h,c] = np.polyval(coeffarray[nl,:,c],ypos)             
                for w in range(self.imagearray.shape[1]):
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
        for i in range(nlevels):
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
        for h in range(od_roi.shape[0]):
            if cont2/cont_lim > 0.05:
                cont2 = 0.0
                print(f'Multichannel correction process: {np.trunc(100*cont/cont_lim)}%')
            for w in range(od_roi.shape[1]):
                #fittet_params = optimize.minimize(self.optimization_func1,
                #                                  x0=np.array([1.0]),
                #                                  args=(od_roi[h, w, :], self.OD0[:])).x
                fittet_params = self.mi_optimizacion(np.array([1.0]), od_roi[h, w, :],self.OD0[:])
                if fittet_params[0] > 1.2:
                    fittet_params[0] = 1.2
                if fittet_params[0] < 0.8:
                    fittet_params[0] = 0.8
                self.alpha[h + y[0], w + x[0]] = fittet_params[0]
                for c in [0, 1, 2]:
                    od = od_roi[h, w, c]*fittet_params[0]
                    self.OD[h + y[0], w + x[0], c] = od
                    cont = cont + 1.0
                    cont2 = cont2 + 1.0
        # Genera una imagen con el mapa de corrección multicanal obtenido
        alpha_image = np.array((255 * (self.alpha-0.8) / 0.4).astype(np.uint8))
        im_alpha = Image.fromarray(alpha_image)
        imname = 'MultiChannelMap_' + self.imagefilename
        im_alpha.save(self.workingdir + imname)

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
        dr = self.Calibration.Get_DoseFromODnet(nod[0], 0)
        dg = self.Calibration.Get_DoseFromODnet(nod[1], 1)
        db = self.Calibration.Get_DoseFromODnet(nod[2], 2)
        f = 3*np.power(dr - dg, 2) + np.power(dr - db, 2) + np.power(dg - db, 2)
        return f

    # Método de corrección multicanal con dos parámetros (alfa, beta).
    # Basado en el método desarrollado por Damian.
    # En desarrollo...
    def multichannel_correction2(self):
        from scipy import optimize
        import ctypes  # An included library with Python install.
        ctypes.windll.user32.MessageBoxW(0, "Select area for multichannel correction", "", 0)
        rs = self.SelectRectangle()
        x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
        y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
        nod_roi = self.OD[y[0]:(y[1] + 1), x[0]:(x[1] + 1), 0:3]
        nod_roi = nod_roi - self.OD0
        self.alpha = np.zeros([self.OD.shape[0], self.OD.shape[1]])
        self.beta = np.zeros([self.OD.shape[0], self.OD.shape[1]])
        lim = 0.2
        cont = 0.0
        cont2 = 0.0
        cont_lim = nod_roi.shape[0] * nod_roi.shape[1] * nod_roi.shape[2]
        for h in range(nod_roi.shape[0]):
            if cont2/cont_lim > 0.05:
                cont2 = 0.0
                print(f'Multichannel correction process: {np.trunc(100*cont/cont_lim)}%')
            for w in range(nod_roi.shape[1]):
                fittet_params = optimize.minimize(self.optimization_func2a,
                                                  x0=np.array([0.01, 0.01]),
                                                  args=nod_roi[h, w, :]).x
                if fittet_params[0] > 0.2:
                    fittet_params[0] = 0.2
                if fittet_params[1] > 0.2:
                    fittet_params[1] = 0.2
                if fittet_params[0] < -0.2:
                    fittet_params[0] = -0.2
                if fittet_params[1] < -0.2:
                    fittet_params[1] = -0.2
                self.alpha[h + y[0], w + x[0]] = fittet_params[0]
                self.beta[h + y[0], w + x[0]] = fittet_params[1]
                for c in [0, 1, 2]:
                    od = nod_roi[h, w, c] + self.OD0[c]
                    od = self.OD0[c] + (od - (1 + fittet_params[1])*self.OD0[c])/(1 + fittet_params[0])
                    self.OD[h + y[0], w + x[0], c] = od
                    cont = cont + 1.0
                    cont2 = cont2 + 1.0
        alpha_image = np.array((255 * self.alpha / 0.4 + 127.5).astype(np.uint8))
        im_alpha = Image.fromarray(alpha_image)
        imname = 'AlphaMap_' + self.imagefilename
        im_alpha.save(self.workingdir + imname)
        beta_image = np.array((255 * self.beta / 0.4 + 127.5).astype(np.uint8))
        im_beta = Image.fromarray(beta_image)
        imname = 'BetaMap_' + self.imagefilename
        im_beta.save(self.workingdir + imname)
        #return self.FScorr


    def optimization_func2a(self, params, odnet):

        nod = (odnet - params[1] * self.OD0)/(1 + params[0])
        dr = self.Calibration.Get_DoseFromODnet(nod[0], 0)
        dg = self.Calibration.Get_DoseFromODnet(nod[1], 1)
        db = self.Calibration.Get_DoseFromODnet(nod[2], 2)
        f = np.power(dr - dg, 2) + np.power(dr - db, 2) + np.power(params[0], 2) + np.power(params[1], 2)
        return f

    def optimization_func2b(self, params, odnet):

        nod = (odnet - params[1] * self.OD0)/(1 + params[0])
        dk = np.array([self.Calibration.Get_DoseFromODnet(odnet[c], c) for c in [0, 1, 2]])
        d = np.array([self.Calibration.Get_DoseFromODnet(nod[c], c) for c in [0, 1, 2]])
        delta_d = np.array([self.Calibration.Get_DerivateFromODnet(nod[c], c) for c in [0, 1, 2]])
        delta_d = delta_d*(params[0]*odnet + params[1]*self.OD0)/(1 + params[0])
        f = 0.0
        unc = np.array([1.0, 1.0, 10.0])
        for c in [0, 1, 2]:
            f = f + np.power(d[c] + delta_d[c] - dk[c], 2)/unc[c]
        #f = f + 0.001*np.power(params[0], 2) + 0.001*np.power(params[1], 2)
        f = f + np.power(d[0]-d[1], 2) + 0.1*np.power(d[0]-d[2], 2)
        return f
        

    