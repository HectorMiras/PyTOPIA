import ntpath
import numpy as np
from lmfit import Model
from pathlib import Path
import FilmDoseClass
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import tkinter as tk
from tkinter import simpledialog
from tkinter import filedialog
from pynverse import inversefunc


class FilmCalibration:
    
    def __init__(self, inputFile=None):
        
        # Parametros iniciales
        self.PercolParam_n = np.array([3, 3, 3])
        self.PercolParam_A = np.array([7, 7, 7])
        self.PercolParam_k = np.array([47, 47, 47])
        self.DevicParam_A = np.array([7.0, 7.0, 7.0])
        self.DevicParam_B = np.array([50.0, 50.0, 50.0])
        self.DevicParam_n = np.array([3.0, 3.0, 3.0])
        self.PV = None
        self.OD0 = None
        self.ODnet = None
        self.PVstd_dev = None
        self.Weights = None
        self.imagefilename = None

        if inputFile is not None:
            self.workingdir = ntpath.dirname(inputFile) + '\\'
            text_extensions = ['.txt', '.dat', '.text']
            tif_extensions = ['.tif', '.tiff']
            if Path(inputFile).suffix in text_extensions:
                # open calibration from text file
                self.calibration_from_text_file(inputFile)
            elif Path(inputFile).suffix in tif_extensions:
                # generate calibration from image
                self.calibration_from_image(inputFile)
                print('Generate Calibration from image')
            else:
                # file with invalid extension
                print("The file has an invalid extension: " + Path(inputFile).suffix)

    def calibration_from_text_file(self, inputFile):
        text_extensions = ['.txt', '.dat', '.text']
        if Path(inputFile).suffix in text_extensions:
            # open calibration from text file
            self.textfilename = ntpath.basename(inputFile)
            try:
                # Fichero que contiene valores de pixel medios y desviaciones estándar
                array_txt = np.loadtxt(inputFile, usecols=(0, 2, 4, 6), skiprows=0)
                self.D = np.loadtxt(inputFile, usecols=(6), skiprows=0)
                if self.D[-1] > 100:
                    self.D[:] = self.D[:]/100
                self.PV = np.loadtxt(inputFile, usecols=(0, 2, 4), skiprows=0)
                self.PVstd_dev = np.loadtxt(inputFile, usecols=(1, 3, 5), skiprows=0)
                self.obtain_ODvalues_from_PV()

            except:
                # Fichero solo contiene valores de pixel medios
                array_txt = np.loadtxt(inputFile,usecols=(0, 1, 2, 3), skiprows=0)
                self.D = np.loadtxt(inputFile,usecols=(3), skiprows=0)
                if self.D[-1] > 100:
                    self.D[:] = self.D[:]/100
                self.PV = np.loadtxt(inputFile,usecols=(0, 1, 2), skiprows=0)
                self.obtain_ODvalues_from_PV()
            self.OptimizeDevicParameters()
        else:
            # file with invalid extension
            print("The file has an invalid extension: " + Path(inputFile).suffix)

    def calibration_from_image(self, inputFile):
        tif_extensions = ['.tif', '.tiff']
        self.workingdir = ntpath.dirname(inputFile) + '\\'
        if Path(inputFile).suffix in tif_extensions:
            # generate calibration from image
            self.imagefilename = ntpath.basename(inputFile)
            self.generate_calibration_from_image()
            print('Generate Calibration from image')
        else:
            # file with invalid extension
            print("The file has an invalid extension: " + Path(inputFile).suffix)

    def calibration_from_arrays(self, pv_array, dose_array, pv_std_dev_array=None):
        self.PV = pv_array
        self.D = dose_array
        if pv_std_dev_array is not None:
            self.PVstd_dev = pv_std_dev_array
        self.obtain_ODvalues_from_PV()
        self.OptimizeDevicParameters()

    def obtain_ODvalues_from_PV(self):
        self.ODnet = np.log10(self.PV[0] / self.PV)
        self.OD0 = np.log10(65535 / self.PV[0])
        if self.PVstd_dev is not None:
            self.Weights = np.power(1.0 / np.log10(self.PV[0] / self.PVstd_dev), 2)
            self.Weights[:, 0] = self.Weights[:, 0] / np.sum(self.Weights[:, 0])
            self.Weights[:, 1] = self.Weights[:, 1] / np.sum(self.Weights[:, 1])
            self.Weights[:, 2] = self.Weights[:, 2] / np.sum(self.Weights[:, 2])
        else:
            self.Weights = np.full(self.PV.shape, 1.0)

    def Get_DoseFromODnet_(self, x, c):
        n = self.PercolParam_n[c]
        A = self.PercolParam_A[c]
        k = self.PercolParam_k[c]
        #return self.PercolCalInv(x, A, k, n)
        return self.DevicCalFunc(x, A, k, n)

    def Get_DoseFromODnet(self, x, c):
        if isinstance(x, (list, np.ndarray)):
            x[x < 0] = 0.0
        else:
            if x < 0:
                x = 0.0
        A = self.DevicParam_A[c]
        B = self.DevicParam_B[c]
        n = self.DevicParam_n[c]
        #return self.PercolCalInv(x, A, k, n)
        return self.DevicCalFunc(x, A, B, n)

    def Get_DerivateFromODnet(self, x, c):
        n = self.PercolParam_n[c]
        A = self.PercolParam_A[c]
        k = self.PercolParam_k[c]
        return self.PercolCalDerivInv(x, A, k, n)


    def GetODnetFromDose_(self, x, c):
        n = self.PercolParam_n[c]
        A = self.PercolParam_A[c]
        k = self.PercolParam_k[c]
        return self.PercolCal(x, A, k, n)

    def GetODnetFromDose(self, x, c):
        n = self.DevicParam_n[c]
        A = self.DevicParam_A[c]
        B = self.DevicParam_B[c]
        return self.DevicInvCalFunc(x, A, B, n)

    def PercolCal(self, x, A, k, n):
        # with np.nditer(x, op_flags=['readwrite']) as it:
        #    for d in it:
        #        if d<0:
        #            d[...] = 0
        y = A * (1 - np.power((1 + x / k), -1 * n))
        return y

    def PercolCalInv(self, x, A, k, n):
        exp = -1 / n
        y = k * (np.power((1 - x / A), exp) - 1)
        return y

    def PercolCalDerivInv(self, x, A, k, n):
        exp = -1 / n
        y = k * (np.power((1 - x / A), (exp - 1)) - 1) / (n * A)
        return y

    def DevicCalFunc(self, x, A, B, n):
        y = A*x + B*np.power(x, n)
        if isinstance(y, (list, np.ndarray)):
            y[y < 0] = 0.0
        else:
            if y < 0:
                y = 0.0
        return y

    def DevicInvCalFunc(self, y, A, B, n):
        f = (lambda x: A*x + B*np.power(x, n))
        od = inversefunc(f, y, domain=[0, 100])
        return od

    def OptimizePercolParameters(self):
        for c in [0, 1, 2]:
            fmodel = Model(self.PercolCal)
            # create parameters -- these are named from the function arguments --
            # giving initial values
            params = fmodel.make_params(A=self.PercolParam_A[c], k=self.PercolParam_k[c], n=self.PercolParam_n[c])
            # fix n:
            #params['n'].vary = False
            # limit parameter values
            params['A'].min = 0
            params['k'].min = 0
            result = fmodel.fit(self.ODnet[:,c], params, x=self.D, weights = self.Weights[:,c])
            self.PercolParam_A[c] = result.params['A']
            self.PercolParam_k[c] = result.params['k']
            self.PercolParam_n[c] = result.params['n']
            print("Fit results for channel {}:".format(c))
            # result.plot_fit()
            print(result.fit_report())
            print('')

    def OptimizeDevicParameters(self):
        for c in [0, 1, 2]:
            fmodel = Model(self.DevicCalFunc)
            # create parameters -- these are named from the function arguments --
            # giving initial values
            params = fmodel.make_params(A=self.DevicParam_A[c], B=self.DevicParam_B[c], n=self.DevicParam_n[c])
            # fix n:
            # params['n'].vary = False
            # limit parameter values
            params['A'].min = 0
            params['A'].max = 100
            params['B'].min = 0
            params['B'].max = 500
            params['n'].min = 1
            params['n'].max = 5
            result = fmodel.fit(self.D, params, x=self.ODnet[:, c])
            self.DevicParam_A[c] = result.params['A']
            self.DevicParam_B[c] = result.params['B']
            self.DevicParam_n[c] = result.params['n']
            print("Fit results for channel {}:".format(c))
            # result.plot_fit()
            print(result.fit_report())
            print('')
        self.show_calibration_plot()
    
    def generate_calibration_from_image(self):
        doseList = []
        roiList = []
        level_count = 0
        filmDose_obj = FilmDoseClass.FilmDose(self.workingdir+self.imagefilename)
        application_window = tk.Tk()
        application_window.withdraw()
        while True:
            answer = simpledialog.askstring("Input", f"Introduce dose (Gy) for level {level_count+1}.",parent=application_window)
            if(answer==None):
                break
            if(np.float(answer)<0 or answer==''): 
                break
            doseList.append(np.float(answer))
            rs = filmDose_obj.SelectRectangle()
            x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
            y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
            ROI_OD = np.array(filmDose_obj.OD[y[0]:(y[1]+1),x[0]:(x[1]+1),0:3])
            roiList.append(ROI_OD)
            level_count = level_count + 1
        self.D = np.zeros(level_count)
        self.OD0 = np.mean(roiList[0],axis=(0,1))
        self.ODnet = np.ones((level_count,3)).astype(np.float)
        self.PV = np.ones((level_count,3)).astype(np.float)
        self.PVstd_dev = np.ones((level_count,3)).astype(np.float)
        OD_std_devs = np.ones((level_count,3)).astype(np.float)
        for i in range(level_count):
            self.D[i] = doseList[i]           
            self.ODnet[i] = np.mean(roiList[i],axis=(0,1))-self.OD0
            OD_std_devs[i] = np.std(roiList[i],axis=(0,1))
            roi = np.array(np.power(10,roiList[i]))
            self.PV[i] = np.mean(65535.0/roi,axis=(0,1))
            self.PVstd_dev[i] = np.std(65535.0/roi,axis=(0,1))
        application_window.destroy()
        self.Weights = np.power(1.0/OD_std_devs,2)
        self.Weights[:,0] = self.Weights[:,0] / np.sum(self.Weights[:,0])
        self.Weights[:,1] = self.Weights[:,1] / np.sum(self.Weights[:,1])
        self.Weights[:,2] = self.Weights[:,2] / np.sum(self.Weights[:,2])
        #self.OptimizePercolParameters()
        self.OptimizeDevicParameters()
        self.show_calibration_plot()

    def save_calibration_text_file(self,calibrationfile):
        file_obj = open(calibrationfile,'w')
        for level in range(len(self.D)):
            L = ''
            for c in [0,1,2]:
                L = L + '{:.2f}'.format(self.PV[level,c]) + ' '
                L = L + '{:.2f}'.format(self.PVstd_dev[level,c]) + ' '
            L = L  + '{:.4f}'.format(self.D[level]) + '\n'
            file_obj.write(L)
        file_obj.close()
    
    def show_calibration_plot(self):
        fig = plt.figure()
        plt.plot(self.ODnet[:,0],self.D, 'ro')
        x = np.arange(0,np.max(self.ODnet[:,0]*1.1),0.01)
        y = self.Get_DoseFromODnet(x,0)
        plt.plot(x,y)
        plt.plot(self.ODnet[:,1],self.D, 'go')
        x = np.arange(0,np.max(self.ODnet[:,1]*1.1),0.01)
        y = self.Get_DoseFromODnet(x,1)
        plt.plot(x,y)
        plt.plot(self.ODnet[:,2],self.D, 'bo')
        x = np.arange(0,np.max(self.ODnet[:,2]*1.1),0.01)
        y = self.Get_DoseFromODnet(x,2)
        plt.plot(x,y)
        plt.show(block=True)



        