#!/usr/bin/env python

### https://physics.nist.gov/PhysRefData/ASD/lines_form.html 

## 
# https://mathematica.stackexchange.com/questions/85990/how-to-plot-an-emission-spectrum
# https://gist.github.com/error454/65d7f392e1acd4a782fc
# https://www.nist.gov/pml/atomic-spectra-database
## 

import os
import sys
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

def wavelength_nm_to_eV(wavelength_nm):
	return 1239.8 / wavelength_nm # eV

def eV_to_wavelength_nm(eV):
	return eV/1239.8 # nm 

def gaussian(x, mean, std_dev):
    return np.exp(-((x - mean) ** 2) / (2 * std_dev ** 2)) / (std_dev * np.sqrt(2 * np.pi))

def wavelength_to_rgb(wavelength, gamma=0.8):
	'''
	copied from https://gist.github.com/error454/65d7f392e1acd4a782fc

	This converts a given wavelength of light to an 
	approximate RGB color value. The wavelength must be given
	in nanometers in the range from 380 nm through 750 nm
	(789 THz through 400 THz).
	Based on code by Dan Bruton
	ttp://www.physics.sfasu.edu/astro/color/spectra.html
	'''
	if 380 <= wavelength <= 440:
		attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
		red = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
		green = 0
		blue = 1
	elif 440 < wavelength <= 490:
		red = 0
		green = (wavelength - 440) / (490 - 440)
		blue = 1
	elif 490 < wavelength <= 510:
		red = 0
		green = 1
		blue = -1 * (wavelength - 510) / (510 - 490)
	elif 510 < wavelength <= 580:
		red = (wavelength - 510) / (580 - 510)
		green = 1
		blue = 0
	elif 580 < wavelength <= 645:
		red = 1
		green = -1 * (wavelength - 645) / (645 - 580)
		blue = 0
	elif 645 < wavelength <= 750:
		attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
		red = (1.0 * attenuation) ** gamma
		green = 0
		blue = 0
	else:
		red = 0
		green = 0
		blue = 0

	return (red, green, blue)


def grating_spectrum(df_lines,
	fname_pdf='spectrum.pdf',fname_jpg='spectrum.jpg',
	linewidthfactor=8.0):
	#xmin=1.65,xmax=3.26):
	#xmin=380,xmax=750):

	fig, ax = plt.subplots(figsize=(10, 1.5))
	#energy_min = 1.65
	#energy_max = 3.26
	energy_min = 1.4
	energy_max = 3.4
	ax.set_xlim(energy_min, energy_max)
	ax.set_ylim(0, 1)
	ax.set_aspect(0.15)
	plt.subplots_adjust(left=0.01,right=0.99,bottom=0.01,top=0.99)
	#plt.xlabel('Energy (eV)')
	#plt.subplots_adjust(left=0.0, right=0.0, bottom=0.0, top=0.0)
	#ax.set_aspect('auto')
	ax.set_facecolor('black')
	ax.spines['bottom'].set_color('white')
	ax.spines['bottom'].set_linewidth(0.5)
	ax.tick_params(axis='x', colors='white')
	ax.xaxis.label.set_color('white')
	ax.yaxis.set_visible(False)

	for index, row in df_lines.iterrows():
		wavelength_nm = float(row['wavelength_nm'])
		energy_eV = wavelength_nm_to_eV(wavelength_nm)
		intensity = linewidthfactor*np.sqrt(float(row['intens_relative']))
		print(wavelength_nm,energy_eV,intensity)
		#intensity = linewidthfactor*np.float(float(row['intens_relative']))
		color = wavelength_to_rgb(wavelength_nm)
		#ax.plot([wavelength, wavelength], [0, 1], color=color, linewidth=1.0*intensity)
		ax.plot([energy_eV, energy_eV], [0, 1], color=color, linewidth=1.0*intensity)
    
	plt.savefig(fname_pdf, format='pdf', facecolor='black')
	plt.savefig(fname_jpg, format='jpg', facecolor='black')	

ef grating_spectrum_gauss(df_lines,
	fname_pdf='spectrum_gauss.pdf',fname_jpg='spectrum_gauss.jpg',
	linewidthfactor=8.0):

	fig, ax = plt.subplots(figsize=(10, 1.5))
	energy_min = 1.4
	energy_max = 3.4
	ax.set_xlim(energy_min, energy_max)
	ax.set_ylim(0, 1)
	ax.set_aspect(0.15)
	plt.subplots_adjust(left=0.01,right=0.99,bottom=0.01,top=0.99)
	#plt.xlabel('Energy (eV)')
	#plt.subplots_adjust(left=0.0, right=0.0, bottom=0.0, top=0.0)
	#ax.set_aspect('auto')
	ax.set_facecolor('black')
	ax.spines['bottom'].set_color('white')
	ax.spines['bottom'].set_linewidth(0.5)
	ax.tick_params(axis='x', colors='white')
	ax.xaxis.label.set_color('white')
	ax.yaxis.set_visible(False)

	for index, row in df_lines.iterrows():
		wavelength_nm = float(row['wavelength_nm'])
		energy_eV = wavelength_nm_to_eV(wavelength_nm)
		intensity = linewidthfactor*np.sqrt(float(row['intens_relative']))
		print(wavelength_nm,energy_eV,intensity)
		#intensity = linewidthfactor*np.float(float(row['intens_relative']))
		color = wavelength_to_rgb(wavelength_nm)
		#ax.plot([wavelength, wavelength], [0, 1], color=color, linewidth=1.0*intensity)
		ax.plot([energy_eV, energy_eV], [0, 1], color=color, linewidth=1.0*intensity)
    
	plt.savefig(fname_pdf, format='pdf', facecolor='black')
	plt.savefig(fname_jpg, format='jpg', facecolor='black')	

def read_csv_to_dflines(csv_file):
	df = pd.read_csv(csv_file)

	print(df.columns)
	if 'obs_wl_vac(nm)' in df.columns:
		df["wavelength_nm"] = df["obs_wl_vac(nm)"]
	elif 'obs_wl_air(nm)' in df.columns:
		df["wavelength_nm"] = df["obs_wl_air(nm)"]		

	df["wavelength_nm"] = df["wavelength_nm"].str.strip('"')
	df["wavelength_nm"] = df["wavelength_nm"].str.strip('=')
	df["wavelength_nm"] = df["wavelength_nm"].str.strip('"')
	df["wavelength_nm"] = pd.to_numeric(df["wavelength_nm"],errors="coerce")

	df["intens"] = df["intens"].str.strip('"')
	df["intens"] = df["intens"].str.strip('=')
	df["intens"] = df["intens"].str.strip('"')
	df["intens"] = pd.to_numeric(df["intens"],errors="coerce")

	df_lines = df[['wavelength_nm','intens']]
	df_lines = df_lines.dropna(how='any')

	max_intensity = df_lines['intens'].max()
	df_lines['intens_relative'] = df_lines['intens']/float(max_intensity)

	return df_lines

args = sys.argv
fname_csv = args[1] # e.g., 'data/Na_I.csv'
df_lines = read_csv_to_dflines(fname_csv) 

fname_pdf = os.path.splitext(os.path.basename(fname_csv))[0] + '_spec.pdf'
fname_jpg = os.path.splitext(os.path.basename(fname_csv))[0] + '_spec.jpg'
grating_spectrum(df_lines,fname_pdf=fname_pdf,fname_jpg=fname_jpg)

fname_pdf = os.path.splitext(os.path.basename(fname_csv))[0] + '_spec_gauss.pdf'
fname_jpg = os.path.splitext(os.path.basename(fname_csv))[0] + '_spec_gauss.jpg'
grating_spectrum_gauss(df_lines,fname_pdf=fname_pdf,fname_jpg=fname_jpg)