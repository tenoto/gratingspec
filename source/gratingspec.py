#!/usr/bin/env python

### https://physics.nist.gov/PhysRefData/ASD/lines_form.html 

## 
# https://mathematica.stackexchange.com/questions/85990/how-to-plot-an-emission-spectrum
# https://gist.github.com/error454/65d7f392e1acd4a782fc
# https://www.nist.gov/pml/atomic-spectra-database
## 

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd

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


def grating_spectrum(df_lines,outpdf='spectrum.pdf'):
	fig, ax = plt.subplots(figsize=(12, 2))
	ax.set_xlim(300, 800)
	ax.set_ylim(0, 1)
	ax.set_aspect('auto')
	ax.set_facecolor('black')
	ax.spines['bottom'].set_color('white')
	ax.spines['bottom'].set_linewidth(0.5)
	ax.tick_params(axis='x', colors='white')
	ax.xaxis.label.set_color('white')
	ax.yaxis.set_visible(False)

	for index, row in df_lines.iterrows():
		wavelength = float(row['wavelength'])
		intensity = float(row['intens_relative'])
		color = wavelength_to_rgb(wavelength)
		ax.plot([wavelength, wavelength], [0, 1], color=color, linewidth=1.0*intensity)
    
	plt.savefig(outpdf, format='pdf', facecolor='black')

def read_csv_to_dflines(csv_file):
	df = pd.read_csv(csv_file)

	print(df.columns)
	if 'obs_wl_vac(nm)' in df.columns:
		df["wavelength"] = df["obs_wl_vac(nm)"]
	elif 'obs_wl_air(nm)' in df.columns:
		df["wavelength"] = df["obs_wl_air(nm)"]		

	df["wavelength"] = df["wavelength"].str.strip('"')
	df["wavelength"] = df["wavelength"].str.strip('=')
	df["wavelength"] = df["wavelength"].str.strip('"')
	df["wavelength"] = pd.to_numeric(df["wavelength"],errors="coerce")

	df["intens"] = df["intens"].str.strip('"')
	df["intens"] = df["intens"].str.strip('=')
	df["intens"] = df["intens"].str.strip('"')
	df["intens"] = pd.to_numeric(df["intens"],errors="coerce")

	df_lines = df[['wavelength','intens']]
	df_lines = df_lines.dropna(how='any')

	max_intensity = df_lines['intens'].max()
	df_lines['intens_relative'] = df_lines['intens']/float(max_intensity)

	return df_lines

args = sys.argv
csvfname = args[1] # e.g., 'data/Na_I.csv'
outfname = os.path.splitext(os.path.basename(csvfname))[0] + '_spec.pdf'
df_lines = read_csv_to_dflines(csvfname) 
grating_spectrum(df_lines,outpdf=outfname)
