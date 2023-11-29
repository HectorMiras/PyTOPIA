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
from scipy.interpolate import UnivariateSpline

class FilmCalibration:

    def __init__(self, inputFile=None):

        # Parametros iniciales
        self.spline_nod_to_D = []
        self.spline_D_to_nod = []
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
                    self.D[:] = self.D[:] / 100
                self.PV = np.loadtxt(inputFile, usecols=(0, 2, 4), skiprows=0)
                self.PVstd_dev = np.loadtxt(inputFile, usecols=(1, 3, 5), skiprows=0)
                self.obtain_ODvalues_from_PV()

            except:
                # Fichero solo contiene valores de pixel medios
                array_txt = np.loadtxt(inputFile, usecols=(0, 1, 2, 3), skiprows=0)
                self.D = np.loadtxt(inputFile, usecols=(3), skiprows=0)
                if self.D[-1] > 100:
                    self.D[:] = self.D[:] / 100
                self.PV = np.loadtxt(inputFile, usecols=(0, 1, 2), skiprows=0)
                self.obtain_ODvalues_from_PV()
            # call calibration construction
            self.create_calibration_nod()

            #self.show_calibration_plot()
            #self.show_inv_calibration_plot()
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
        # call calibration construction
        self.create_calibration_nod()

    def generate_calibration_from_image(self):
        doseList = []
        roiList = []
        level_count = 0
        filmDose_obj = FilmDoseClass.FilmDose(self.workingdir + self.imagefilename)

        while True:
            application_window = tk.Tk()
            application_window.withdraw()
            answer = simpledialog.askstring("Input", f"Introduce dose (Gy) for level {level_count + 1}.",
                                            parent=application_window)
            if (answer == None):
                application_window.destroy()
                break
            if (np.float(answer) < 0 or answer == ''):
                application_window.destroy()
                break
            application_window.destroy()
            doseList.append(np.float(answer))
            rs = filmDose_obj.SelectRectangle()
            x = [np.int(rs.corners[0][0]), np.int(rs.corners[0][2])]
            y = [np.int(rs.corners[1][0]), np.int(rs.corners[1][2])]
            ROI_OD = np.array(filmDose_obj.OD[y[0]:(y[1] + 1), x[0]:(x[1] + 1), 0:3])
            roiList.append(ROI_OD)
            level_count = level_count + 1
        self.D = np.zeros(level_count)
        self.OD0 = np.mean(roiList[0], axis=(0, 1))
        self.ODnet = np.ones((level_count, 3)).astype(np.float)
        self.PV = np.ones((level_count, 3)).astype(np.float)
        self.PVstd_dev = np.ones((level_count, 3)).astype(np.float)
        OD_std_devs = np.ones((level_count, 3)).astype(np.float)
        for i in range(level_count):
            self.D[i] = doseList[i]
            self.ODnet[i] = np.mean(roiList[i], axis=(0, 1)) - self.OD0
            OD_std_devs[i] = np.std(roiList[i], axis=(0, 1))
            roi = np.array(np.power(10, roiList[i]))
            self.PV[i] = np.mean(65535.0 / roi, axis=(0, 1))
            self.PVstd_dev[i] = np.std(65535.0 / roi, axis=(0, 1))

        self.Weights = np.power(1.0 / OD_std_devs, 2)
        self.Weights[:, 0] = self.Weights[:, 0] / np.sum(self.Weights[:, 0])
        self.Weights[:, 1] = self.Weights[:, 1] / np.sum(self.Weights[:, 1])
        self.Weights[:, 2] = self.Weights[:, 2] / np.sum(self.Weights[:, 2])

        # Call calibrationconstruction
        self.create_calibration_nod()

        self.show_calibration_plot()

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

    def create_calibration_nod(self):
        for c in range(3):
            self.spline_nod_to_D.append(UnivariateSpline(self.ODnet[:, c], self.D, s=0))
            self.spline_D_to_nod.append(UnivariateSpline(self.D, self.ODnet[:, c], s=0))

    def get_dose_from_nod(self, x, ch=None):
        dim = x.ndim
        x = np.clip(x, 0, None)
        dose_arr = np.array(0.0)
        if dim == 3:
            if ch is None:
                dose_arr = np.array([self.spline_nod_to_D[c](x[:, :, c]) for c in np.arange(3)])
                dose_arr = dose_arr.swapaxes(0, 1).swapaxes(1, 2)
            else:
                dose_arr = self.spline_nod_to_D[ch](x)
        else:
            if ch is None:
                dose_arr = np.array([self.spline_nod_to_D[c](x[c]) for c in np.arange(3)])
            else:
                dose_arr = self.spline_nod_to_D[ch](x)

        return dose_arr

    def get_nod_from_dose(self, x, ch=None):
        dim = x.ndim
        x = np.clip(x, 0, None)
        nod = np.array(0.0)
        if ch is None:
            nod = np.array([self.spline_D_to_nod[c](x[c]) for c in np.arange(3)])
            if nod.ndim == 3:
                nod = nod.swapaxes(0, 1).swapaxes(1, 2)
        else:
            nod = self.spline_D_to_nod[ch](x)

        return nod

    def save_calibration_text_file(self, calibrationfile):
        file_obj = open(calibrationfile, 'w')
        for level in range(len(self.D)):
            L = ''
            for c in [0, 1, 2]:
                L = L + '{:.2f}'.format(self.PV[level, c]) + ' '
                L = L + '{:.2f}'.format(self.PVstd_dev[level, c]) + ' '
            L = L + '{:.4f}'.format(self.D[level]) + '\n'
            file_obj.write(L)
        file_obj.close()

    def show_calibration_plot(self):
        fig = plt.figure()
        plt.plot(self.ODnet[:, 0], self.D, 'ro')
        x = np.arange(0, np.max(self.ODnet[:, 0] * 1.1), 0.01)
        y = self.get_dose_from_nod(x, 0)
        plt.plot(x, y)
        plt.plot(self.ODnet[:, 1], self.D, 'go')
        x = np.arange(0, np.max(self.ODnet[:, 1] * 1.1), 0.01)
        y = self.get_dose_from_nod(x, 1)
        plt.plot(x, y)
        plt.plot(self.ODnet[:, 2], self.D, 'bo')
        x = np.arange(0, np.max(self.ODnet[:, 2] * 1.1), 0.01)
        y = self.get_dose_from_nod(x, 2)
        plt.plot(x, y)
        plt.show(block=True)

    def show_inv_calibration_plot(self):
        fig = plt.figure()
        plt.plot(self.D, self.ODnet[:, 0], 'ro')
        x = np.arange(0, np.max(self.D * 1.1), 0.01)
        y = self.get_nod_from_dose(x, 0)
        plt.plot(x, y)
        plt.plot(self.D, self.ODnet[:, 1], 'go')
        y = self.get_nod_from_dose(x, 1)
        plt.plot(x, y)
        plt.plot(self.D, self.ODnet[:, 2], 'bo')
        y = self.get_nod_from_dose(x, 2)
        plt.plot(x, y)
        plt.show(block=True)