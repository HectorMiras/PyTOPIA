import numpy as np
import FilmDoseClass
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog


def main():
    
    # carga la imagen a procesar
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title='Open Tif Image')
    filmDoseObj = FilmDoseClass.FilmDose(file_path)
    root.destroy()
    
    # Carga los ficheros de flatscanner (1D)
    # selecciona cancelar cuando ya no se quieran cargar más ficheros
    filmDoseObj.LoadFlatScanCorrFiles()
    # Aplica corrección flatscan (1D)
    FSmatrix = filmDoseObj.ApplyFlatScanCorr()

    # Carga calibración desde fichero
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title='Open Calibration file')
    root.destroy()

    filmDoseObj.AddCalibrationFromFile(file_path)
    # Aplica corrección de la calibración
    filmDoseObj.CalibrationCorrection()
    print('Calibration correction factors:')
    print(filmDoseObj.CalibrationCorrectionFactors)

    # Aplica método de corrección multicanal de 1 parámetro
    filmDoseObj.multichannel_correction1()

    # Salva la imagen en dosis.
    file_path = filedialog.asksaveasfilename(title='Save Corrected Dose Image')
    # Pide la dosis máxima a la que corresponderá el valor de pixel máximo
    application_window = tk.Tk()
    application_window.withdraw()
    answer = simpledialog.askstring("Input", "Introduce max dose (Gy) for linearization",
                                    parent=application_window)
    filmDoseObj.SaveDoseImage(file_path, np.float(answer))
    application_window.destroy()

    # Salva imagen en valor de pixel pero con todas las correcciones aplicadas
    #file_path = filedialog.asksaveasfilename(title='Save Corrected Pixel Value Image')
    #filmDoseObj.SaveCorrectedPVImage(file_path)
    
    print('End of program.')

if __name__ == '__main__':
    
    main()
    