import FilmDoseClass
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
import config
import glob
import os


def get_files(directory_path, pattern):
    # Construct the full pattern
    full_pattern = os.path.join(directory_path, pattern)

    # List all files in the directory that match the pattern
    matching_files = glob.glob(full_pattern)

    # return the list of matching files
    files_list = []
    for file in matching_files:
        files_list.append(file)

    return files_list

def main():
    files_list = []

    if config.file_pattern and config.working_dir:
        files_list = get_files(config.working_dir, config.file_pattern)
    else:
        # carga la imagen a procesar
        root = tk.Tk()
        root.withdraw()
        files_list.append(filedialog.askopenfilename(title='Open Tif Image'))
        root.destroy()

    for file_path in files_list:
        print(f'Processing file {file_path}:')
        filmDoseObj = FilmDoseClass.FilmDose(file_path)

        # Carga los ficheros de flatscanner (1D)
        # selecciona cancelar cuando ya no se quieran cargar más ficheros
        if config.ApplyFS:
            filmDoseObj.LoadFlatScanCorrFiles(config.FS_files_path)
            # Aplica corrección flatscan (1D)
            FSmatrix = filmDoseObj.ApplyFlatScanCorr()

        # Carga calibración desde fichero
        if config.ApplyCalibration:
            if config.CalibrationFile:
                cal_file = config.CalibrationFile
            else:
                root = tk.Tk()
                root.withdraw()
                cal_file = filedialog.askopenfilename(title='Open Calibration file')
                root.destroy()

            filmDoseObj.AddCalibrationFromFile(cal_file)

        # Aplica corrección de la calibración
        if config.ApplyCalCorrection:
            filmDoseObj.CalibrationCorrection()
            print('Calibration correction factors:')
            print(filmDoseObj.CalibrationCorrectionFactors)

        # Aplica método de corrección multicanal
        if config.ApplyMultichannel:
            filmDoseObj.multichannel_correction1()

        if config.ApplyEnergyDependence:
            filmDoseObj.EnergyDependence_correction(config.EnDep_kint, config.EnDep_fad)

        # Salva la imagen en dosis.
        if config.SaveDoseFile:
            if config.OutputMaxDose is None:
                application_window = tk.Tk()
                application_window.withdraw()
                # Pide la dosis máxima a la que corresponderá el valor de pixel máximo
                maxdose = simpledialog.askstring("Input", "Introduce max dose (Gy) for linearization",
                                                parent=application_window)
                application_window.destroy()
            else:
                maxdose = config.OutputMaxDose

            # Guarda la imagen

            if config.AutomaticOutput:
                root, extension = os.path.splitext(file_path)
                output_name = os.path.basename(file_path).split('.')[-2]+"-"
                if config.ApplyCalibration:
                    output_name += 'D'
                if config.ApplyCalCorrection:
                    output_name += 'K'
                if config.ApplyFS:
                    output_name += 'S'
                if config.ApplyEnergyDependence:
                    output_name += 'E'
                if config.ApplyMultichannel:
                    output_name += 'M'
                output_name += extension
                output_name = os.path.join(config.working_dir, output_name)
            else:
                application_window = tk.Tk()
                application_window.withdraw()
                output_name = filedialog.asksaveasfilename(title='Save Corrected Dose Image')
                application_window.withdraw()

            filmDoseObj.SaveDoseImage(output_name, float(maxdose))

        if config.SavePVFile:
            if config.AutomaticOutput:
                root, extension = os.path.splitext(file_path)
                output_name = os.path.basename(file_path).split('.')[-2]+"-"
                if config.ApplyCalibration:
                    output_name += 'D'
                if config.ApplyCalCorrection:
                    output_name += 'K'
                if config.ApplyFS:
                    output_name += 'S'
                if config.ApplyEnergyDependence:
                    output_name += 'E'
                if config.ApplyMultichannel:
                    output_name += 'M'
                output_name += "_PV"
                output_name += extension
                output_name = os.path.join(config.working_dir, output_name)
            else:
                application_window = tk.Tk()
                application_window.withdraw()
                output_name = filedialog.asksaveasfilename(title='Save Corrected Pixel Value Image')
                application_window.withdraw()
            filmDoseObj.SaveCorrectedPVImage(output_name)

        print('')
    print('End of program.')


if __name__ == '__main__':
    main()