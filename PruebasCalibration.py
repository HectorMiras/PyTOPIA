import FilmCalibrationClass
import tkinter as tk
from tkinter import filedialog

def main():
    tk.Tk().withdraw()
    file_path = filedialog.askopenfilename(title='Open Calibration (image or text)')
    if(file_path != ''):
        filmCal_obj = FilmCalibrationClass.FilmCalibration(file_path)
        file_path = filedialog.asksaveasfilename(title='Save calibration file')
        filmCal_obj.save_calibration_text_file(file_path)
    print('End')


if __name__ == '__main__':
    
    main()
    