#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file_name.py
# Python 3.8.3
"""
Author: Joseph Bateman
Created: Fri Nov  5 16:01:25 2021
Modified: Fri Nov  5 16:01:25 2021

Description 
-------------

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from astropy.modeling import models, fitting


def get_rgb_d_from_file(filename):
    input_image = Image.open(filename).convert("RGB")
    input_imarray = np.array(input_image)

    X_d = np.array(input_imarray / 65535)

    X_d_transpose = np.transpose(X_d)

    rX_d, gX_d, bX_d = X_d_transpose
    return rX_d, gX_d, bX_d


def get_dose(a, b, c, X_d):
    dose = np.fliplr(np.rot90((a - c * X_d) / (X_d - b), 3))
    return dose


# fit 2d gaussian to get beam size
def fit_2D_gaussian_to_dose(dose_profile):

    y, x = (np.mgrid[:np.shape(dose_profile)[0], :np.shape(dose_profile)[1]])
    g = models.Gaussian2D()
    fitter = fitting.LevMarLSQFitter()
    p = fitter(g, x, y, dose_profile)
    plt.figure(figsize=(8, 2.5))
    plt.subplot(1, 2, 1)
    im1=plt.imshow(dose_profile, cmap="jet", interpolation="nearest", vmin=0, vmax=10)
    plt.title("Raw Dose Map")
    plt.colorbar(im1, label= 'Dose(Gy)')
    plt.xlabel('X pos (pixels)')
    plt.ylabel('Y pos (pixels)')
    plt.subplot(1, 2, 2)
    im2=plt.imshow(p(x, y), cmap="jet", interpolation="nearest", vmin=0, vmax=10)
    plt.title("Gauss Fit Dose Map")
    plt.colorbar(im2, label='Dose (Gy)')
    plt.xlabel('X pos (pixels)')
    plt.ylabel('Y pos (pixels)')
    #plt.subplot(1, 3, 3)
    #plt.imshow(dose_profile - p(x, y), cmap="jet", interpolation="nearest", vmin=-2, vmax=2)
    #plt.title("Residual")
    #plt.colorbar()
    return p


def make_fit_params_table(red_fit, green_fit, blue_fit):
    fit_list = [red_fit, green_fit, blue_fit]
    colours = ["Red", "Green", "Blue"]
    amplitudes = []
    x_means = []
    y_means = []
    x_stddevs = []
    y_stddevs = []
    for i in fit_list:
        amplitudes.append(i.amplitude)
        x_means.append(i.x_mean) 
        y_means.append(i.y_mean) 
        x_stddevs.append(i.x_stddev*25.5/300)
        y_stddevs.append(i.y_stddev*25.4/300)

    fitting_data = pd.DataFrame(
        {
            "Colour": colours,
            "Amplitude": amplitudes,
            "x_mean": x_means,
            "y_mean": y_means,
            "x_stddev": x_stddevs,
            "y_stddev": y_stddevs,
        }
    )
    return fitting_data

a_r = 1.08979454e-02
b_r = 4.31528862e-04
c_r = 4.52890034e+00
  
a_g = 1.30365440e-02
b_g = 2.18523305e-04
c_g = 5.47818342e+00
  
a_b = 1.83975349e-02
b_b = 4.03633386e-04
c_b = 1.10176866e+01

    
rX_d, gX_d, bX_d = get_rgb_d_from_file("/Users/josephbateman/CLEAR/Film/UVic/SFRT_09_11_22/70MeV/70Mev001.tif")

dose_red = get_dose(a_r, b_r, c_r, rX_d)
dose_green = get_dose(a_g, b_g, c_g, gX_d)
dose_blue = get_dose(a_b, b_b, c_b, bX_d)

fitted_red = fit_2D_gaussian_to_dose(dose_red)
fitted_green = fit_2D_gaussian_to_dose(dose_green)
fitted_blue = fit_2D_gaussian_to_dose(dose_blue)

fitting_data_table = make_fit_params_table(fitted_red, fitted_green, fitted_blue)
print(fitting_data_table)

#x mean position green channel
xmean=fitting_data_table.iloc[1,2]*25.4/300


#y mean position green channel
ymean=fitting_data_table.iloc[1,3]*25.4/300

dim = 10

max_dose=np.round(np.mean(dose_green[int(fitting_data_table.iloc[1,3]-dim/2):int(fitting_data_table.iloc[1,3]+dim/2),int(fitting_data_table.iloc[1,2]-dim/2):int(fitting_data_table.iloc[1,2]+dim/2)]),2)

print(max_dose, fitting_data_table.iloc[1,4], fitting_data_table.iloc[1,5])





