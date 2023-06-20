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
import os
import matplotlib.pyplot as plt
from PIL import Image
from astropy.modeling import models, fitting
from scipy.signal import savgol_filter

# class holding fitting and dose measurement methods
class film_analyser:
    # initialise by defining calibration constants
    # and retrieving dose in seperate channels from given file
    def __init__(self, filename):
        # seperate film into rgb channels and return data values
        def get_rgb_d_from_file(filename):
            # open file and convert to Image object in RGB format for analysis
            input_image = Image.open(filename).convert("RGB")
            # convert to numpy array for ease
            input_imarray = np.array(input_image)
            # X_d reweighted because reasons
            X_d = np.array(input_imarray / 65535)
            # transpose carried out to retrieve channels from matrix
            X_d_transpose = np.transpose(X_d)
            # values in r,g,b colour channels retrieved
            rX_d, gX_d, bX_d = X_d_transpose
            return rX_d, gX_d, bX_d

        # method to convert image matrix to dose with calibration constants
        def get_dose(a, b, c, X_d):
            # conversion formula and matrix manipulation to retrieve dose
            dose = np.fliplr(np.rot90((a - c * X_d) / (X_d - b), 3))
            return dose
        
        #ebt3 calibration params up to 25 Gy

        #a_r = 9.08646370e-03
        #b_r = 4.83772917e-04
        #c_r = 3.35869684e+00
        
        #a_g = 1.25584443e-02
        #b_g = 2.13539882e-04
        #c_g = 4.98839944e+00
        
        #a_b = 1.84179296e-02
        #b_b = 3.86486899e-04
        #c_b = 1.10265417e+01

        #ebt3 full calibration params

       # a_r = 1.08054260e-02
       # b_r = 3.96464552e-04
       # c_r = 4.03747782e+00

       # a_g = 1.53583642e-02
       # b_g = 9.04451292e-05
       # c_g = 6.24521113e+00

       # a_b = 2.25123199e-02
       # b_b = 2.75485087e-04
       # c_b = 1.36316506e+01
       
       #mdv3 new calibration params

        #a_r = 7.21071436e-02
        #b_r = 4.67275120e-04
        #c_r = 2.38511751e+01

        #a_g = 1.11306451e-01
        #b_g = 3.18670259e-04
        #c_g = 3.79674212e+01

        #a_b = 2.32990669e-01
        #b_b = 7.33015188e-04
        #c_b = 9.54995622e+01
        
        #mdv3 old calibration params

        #a_r = 7.33992768e-02
        #b_r = 4.73098274e-04
        #c_r = 2.47012235e+01

        #a_g = 1.16746035e-01
        #b_g = 2.36675318e-04
        #c_g = 4.08090283e+01

        #a_b = 2.77850507e-01
        #b_b = 4.96471261e-04
        #c_b = 1.17878125e+02
        
        #EBT3 New calib
        
        a_r = 1.08979454e-02
        b_r = 4.31528862e-04
        c_r = 4.52890034e+00
          
        a_g = 1.30365440e-02
        b_g = 2.18523305e-04
        c_g = 5.47818342e+00
          
        a_b = 1.83975349e-02
        b_b = 4.03633386e-04
        c_b = 1.10176866e+01

        
        # retrieve red, green and blue image data from file
        rX_d, gX_d, bX_d = get_rgb_d_from_file(filename)
        # retrieve dose in each colour channel
        self.dose_red = get_dose(a_r, b_r, c_r, rX_d)
        self.dose_green = get_dose(a_g, b_g, c_g, gX_d)
        self.dose_blue = get_dose(a_b, b_b, c_b, bX_d)


# path must end with a forward slash
# function to loop through and analyse every film image in a directory
def analyse_all_films(path_to_directory_holding_images):
    # call film analyser object to retrieve dose and get gaussian fitting data
    def get_dose_green(filename):
        # initialise film analyser object
        film = film_analyser(filename)
        # get doses in colour channels
        dose_red, dose_green, dose_blue = film.dose_red, film.dose_green, film.dose_blue
   
        # return  green dose for further processing
        return dose_green

    # return max dose and save to file
    def get_max_dose(dose_green):

     
        dim = 10
        film_y_size, film_x_size = np.shape(dose_green)
      
        
        # find max dose
        max_dose = np.round(
            np.mean(
                dose_green[
                    int(227) : int(
                        229
                    ),
                    int(241) : int(
                        243
                    ),
                ]
            ),
            2,
        )
        # print max dose
        print(max_dose)
        # save results to .csv file
        with open("/Users/josephbateman/CLEAR/Film/UVic/SFRT_09_11_22/70MeVcodecheck.csv", "a", newline="") as csv_file:
            print(
                max_dose,
                filename,
                file=csv_file,
                sep=",",
            )

    # get directory in bytes format
    directory = os.fsencode(path_to_directory_holding_images)
    # begin loop for analysing every file in folder
    for file in os.listdir(directory):

        # decode from bytes into string
        file_str = file.decode("utf-8")
        if file_str[0] != ".":
            filename = file_str
            # add path to directory to filename to point correctly
            file_str = path_to_directory_holding_images + file_str
            # generate and retrieve fitting table and green_dose
            dose_green = get_dose_green(
                file_str
            )
            print("File: " + filename)
            # retrieve and print max_dose
            get_max_dose(dose_green)


# main function
def main():
    # carry out full analysis of films in chosen directory
    analyse_all_films("/Users/josephbateman/CLEAR/Film/UVic/SFRT_09_11_22/70MeV/")


main()
