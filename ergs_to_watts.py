import pandas as pd
import astropy.units as u
import astropy.constants as C
import numpy as np

# Code to convert data from Janskys to u.erg/u.s/u.cm**2/u.AA to then have Brewster G935H convert to to(u.W/u.m**2/u.um)

#Read in the data
path_to_spectra= '/Users/eileengonzales/Box/Research/JWST_ross458c/Spectra/'
df= pd.read_csv(path_to_spectra+"Wolf1130C.txt",sep="\s+", comment="#", names=['w','f','e'])

#convert Janskys to Ergs/s/cm^2/um. Need to add a factor of 1/lambda^2 to account for nu to lambda.
Flux_new= ((df['f'].values*(u.Jy)* C.c/(df['w'].values*u.micron)**2)).to(u.erg/u.s/u.cm**2/u.AA).value
flux_err_new =((df['e'].values*(u.Jy)* C.c/(df['w'].values*u.micron)**2)).to(u.erg/u.s/u.cm**2/u.AA).value

#Create a new dataframe with the converted values and write out to a csv file.
df_new= pd.DataFrame()
df_new['w'] = df['w']
df_new['f'] = Flux_new
df_new['e'] = flux_err_new

#Write to CSV file without the header. not sure how to comment it out instead.
df_new.to_csv(path_to_spectra+"converted_units_Wolf1130C.txt",index=False,sep='\t', header=None)
#FYI conversion between to(u.erg/u.s/u.cm**2/u.AA) to W/m**2/um is just to mulitply by 10! (see G395H template for example)









