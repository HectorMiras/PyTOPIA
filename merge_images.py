import numpy as np
from PIL import Image
import os
import tifffile #https://pypi.org/project/tifffile/

# Image pattern
im_pattern = "Film25_Calibracion*"
wdir = r"C:\datos\DosimetriaPelicula\Imagenes\FLASH_24-11-2023"

# List of image file paths
image_files = ['Film24_CalibracionFotones001.tif',
               'Film24_CalibracionFotones002.tif',
               'Film24_CalibracionFotones003.tif',
               'Film24_CalibracionFotones004.tif',
               'Film24_CalibracionFotones005.tif']

# Initialize an array to hold the sum of all images
sum_images = None

# Iterate over image files to calculate the sum
for image_file in image_files:
    # Open image and convert to numpy array
    img = np.array(tifffile.imread(os.path.join(wdir, image_file)))
    tif = tifffile.TiffFile(os.path.join(wdir, image_file))
    tag = tif.pages[0].tags['XResolution']
    dpi = tag.value[0] / tag.value[1]
    # If it's the first image, initialize the sum_images array
    if sum_images is None:
        sum_images = np.array(img, dtype=np.float64)
    else:
        # Add current image to the sum_images array
        sum_images += np.array(img, dtype=np.float64)

# Calculate the average
average_image = sum_images / len(image_files)

# Convert back to uint16
average_image = np.clip(average_image, 0, 65535)  # Ensure values are in 16-bit range
average_image = average_image.astype(np.uint16)

# Save the result
tifffile.imwrite(os.path.join(wdir,'average_image.tif'), average_image, resolution=(dpi, dpi))
