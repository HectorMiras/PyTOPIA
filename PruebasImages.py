import numpy as np
import FilmDoseClass
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog


def main():

    ApplyFS = True
    ApplyCalibration = True
    ApplyCalCorrection = True
    ApplyMultichannel = True
    GuardaDosis = True
    GuardaPV = False

    # carga la imagen a procesar
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title='Open Tif Image')
    filmDoseObj = FilmDoseClass.FilmDose(file_path)
    root.destroy()
    
    # Carga los ficheros de flatscanner (1D)
    # selecciona cancelar cuando ya no se quieran cargar más ficheros
    if ApplyFS:
        filmDoseObj.LoadFlatScanCorrFiles()
        # Aplica corrección flatscan (1D)
        FSmatrix = filmDoseObj.ApplyFlatScanCorr()

    # Carga calibración desde fichero
    if ApplyCalibration:
        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename(title='Open Calibration file')
        root.destroy()
        filmDoseObj.AddCalibrationFromFile(file_path)
        dosis = filmDoseObj.get_dose()
    # Aplica corrección de la calibración
    if ApplyCalCorrection:
        filmDoseObj.CalibrationCorrection()
        print('Calibration correction factors:')
        print(filmDoseObj.CalibrationCorrectionFactors)

    # Asigna parametros alpha y beta y sus sigmas
    #filmDoseObj.Calibration.AlphaCal = 0.0
    #filmDoseObj.Calibration.SigmaAlphaCal = 1.0/np.sqrt(1.3083e4)
    #filmDoseObj.Calibration.BetaCal = 0.0
    #filmDoseObj.Calibration.SigmaBetaCal = 1.0/np.sqrt(7.5337e3)

    #filmDoseObj.Calibration.DevicParam_A = np.array([6.702, 11.3167, 28.6477])
    #filmDoseObj.Calibration.DevicParam_B = np.array([47.5963, 43.9465, 160.7376])
    #filmDoseObj.Calibration.DevicParam_n = np.array([3.04, 2.47, 2.56])

    # Aplica método de corrección multicanal de Damian
    if ApplyMultichannel:
        filmDoseObj.multichannel_correction_Dam()

    # Salva la imagen en dosis.
    if GuardaDosis:
        application_window = tk.Tk()
        application_window.withdraw()
        # Pide la dosis máxima a la que corresponderá el valor de pixel máximo
        answer = simpledialog.askstring("Input", "Introduce max dose (Gy) for linearization",
                                        parent=application_window)
        application_window.destroy()
        # Guarda la imagen
        application_window = tk.Tk()
        application_window.withdraw()
        file_path = filedialog.asksaveasfilename(title='Save Corrected Dose Image')
        application_window.withdraw()
        filmDoseObj.SaveDoseImage(file_path, np.float(answer))

    if GuardaPV:
        application_window = tk.Tk()
        application_window.withdraw()
        file_path = filedialog.asksaveasfilename(title='Save Corrected Pixel Value Image')
        application_window.withdraw()
        filmDoseObj.SaveCorrectedPVImage(file_path)
    
    print('End of program.')

if __name__ == '__main__':
    
    main()
    