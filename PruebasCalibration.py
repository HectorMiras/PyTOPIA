import FilmCalibrationClass_Dam
import tkinter as tk
from tkinter import filedialog

def main():
    appwindow = tk.Tk()
    appwindow.withdraw()
    file_path = filedialog.askopenfilename(title='Open Calibration (image or text)')
    appwindow.destroy()
    if(file_path != ''):
        filmCal_obj = FilmCalibrationClass_Dam.FilmCalibration(file_path)
        if(filmCal_obj.imagefilename is not None):
            appwindow = tk.Tk()
            appwindow.withdraw()
            file_path = filedialog.asksaveasfilename(title='Save calibration file')
            if len(file_path) > 0:
                filmCal_obj.save_calibration_text_file(file_path)
            appwindow.destroy()
    print('End')


if __name__ == '__main__':
    
    main()
    