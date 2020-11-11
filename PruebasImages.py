import numpy as np
import FilmDoseClass
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog


def main():

    ApplyFS = False
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
        # Asigna parametros alpha y beta y sus sigmas
        # filmDoseObj.Calibration.AlphaCal = 5.74357e-5
        #filmDoseObj.Calibration.AlphaCal = 0
        #filmDoseObj.Calibration.SigmaAlphaCal = 0.00901
        # filmDoseObj.Calibration.BetaCal = -0.000115175
        #filmDoseObj.Calibration.BetaCal = 0
        #filmDoseObj.Calibration.SigmaBetaCal = 0.0113175

        #filmDoseObj.Calibration.DevicParam_A = np.array([7.02691762, 11.04462102, 29.33266638])
        #filmDoseObj.Calibration.DevicParam_B = np.array([47.52669613, 43.42677272, 170.73321208])
        #filmDoseObj.Calibration.DevicParam_n = np.array([3.1027201, 2.43657698, 2.64754729])
        #filmDoseObj.Calibration.Sigma_A = np.array([0.15806952, 0.09939978, 0.41270207])
        #filmDoseObj.Calibration.Sigma_B = np.array([0.38155465, 0.21216226, 2.47785548])

    # Aplica corrección de la calibración
    if ApplyCalCorrection:
        filmDoseObj.CalibrationCorrection()
        print('Calibration correction factors:')
        print(filmDoseObj.CalibrationCorrectionFactors)



    # Aplica método de corrección multicanal de Damian
    if ApplyMultichannel:
        filmDoseObj.multichannel_correction_Dam()
        #filmDoseObj.multichannel_correction1()

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
    