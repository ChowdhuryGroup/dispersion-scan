# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog


def filedlg(title=''):
    """
    Pulls up a dialog box to select multiple files. 'title' argument only labels the dialog box
    """
   
    root=tk.Tk()
    files=list(filedialog.askopenfilenames(title='MANDORLA - '+title))
    root.destroy()
    return files

def read_pgm(filename, byteorder='>'):
    """Return image data from a raw PGM file as numpy array.
    Format specification: http://netpbm.sourceforge.net/doc/pgm.html
    """
    with open(filename, 'rb') as f:
        buffer = f.read()
    try:
        header, width, height, maxval = re.search(
            b"(^P5\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n]\s)*)", buffer).groups()
    except AttributeError:
        raise ValueError("Not a raw PGM file: '%s'" % filename)
    return np.frombuffer(buffer,
                            dtype='u1' if int(maxval) < 256 else byteorder+'u2',
                            count=int(width)*int(height),
                            offset=len(header)
                            ).reshape((int(height), int(width)))

def openim(fname,flip_x=False,flip_y=False):
    """"
    A wrapper to open an image and flip the x and/or y axis as needed. Image
    """
    if fname[-3:]=='pgm':
        image=read_pgm(fname)
    else:
        image=plt.imread(fname)
   
    if flip_x:
        image=np.flip(image,axis=1)
    else:
        pass
   
    if flip_y:
        image=np.flip(image,axis=0)
    else:
        pass
   
    plt.imshow(image)
    return image

def fftfilter(image,x_width,y_width):
    """
    Applies a rectangular filter to the Fourier transform passing through low frequencies at the center.
    """
   
    fimage=np.fft.fft2(image)
    filterim=fimage.copy()
    a=x_width
    b=y_width
    filterim[b:-b,a:-a]=0
    newimage=np.abs(np.fft.ifft2(filterim))
   
    return newimage

def findheader(file,delim,maxrows,tries):
    """
    Determines the header for a given file. A brute force method based on whether importing the file throws a ValueError in Numpy.loadtxt
    """
   
    head=0
    headerfound=False
    while head<tries:
        try:
            np.loadtxt(file,delimiter=delim,skiprows=head,max_rows=maxrows)
        except ValueError:
            head +=1
        else:
            headerfound=True
            break
    if headerfound==False:
        print('Header not found after '+str(tries)+' lines.')
    else:
       print('Header was found')
    return head

def findfooter(file,delim,header,maxrows,tries):
    """
    Determines the footer for a given file. A brute force method based on whether importing the file throws a ValueError in Numpy.loadtxt
    """
   
    foot=0
    footerfound=False
    while foot<tries:
        try:
            np.loadtxt(file,delimiter=delim,skiprows=header,max_rows=maxrows-foot)
        except ValueError:
            foot +=1
        else:
            footerfound=True
            break
    if footerfound==False:
        print('Footer not found after '+str(tries)+' lines.')
    else:
       print('Footer was found')
    return foot

def flength(fname):
    """
    Opens a file to determine the number of rows
    """
   
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def finddata(file,delim,tries):
    """
    Determines file size, header, and footer of generic data file. Assumes footer is shorter than header.
    """
   
    totalrows=flength(file)
    header=findheader(file,delim,totalrows-tries,tries)
    maxrows=totalrows-header
    footer=findfooter(file,delim,header,maxrows,tries)
    return maxrows, header, footer

def loadfile(file,delim,tries,row_wise=False):
    """
    Loads data from a file into a numpy array. This method automatically selects out the header and footer with the 'findheader' and 'findfooter' methods
    """
    mr,h,f=finddata(file,delim,tries)
    data=np.loadtxt(file,delimiter=delim,skiprows=h,max_rows=mr-f)
    if row_wise:
        data=data.T
    # loads the data in rows instead of columns instead of columns for easier manipulation
    else:
        pass    
    return data

def opendscan():
    x = filedlg()
    data = loadfile(x[0],"\t", 1)
    spectra = data[0,1:]
    stage = data[2:,0]
    signal = data[2:, 1:]
    Spectra, Stage = np.meshgrid(spectra, stage)
    plot = plt.pcolormesh(Spectra, Stage, signal)

    
    