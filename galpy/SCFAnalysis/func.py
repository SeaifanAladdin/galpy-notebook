import numpy as nu

def meanSquares(raw_dec, fit_dec):
    return (nu.sum((raw_dec - fit_dec)**2)/ len(raw_dec))**.5

def ra_decFit(ra, dec, minRa, maxRa, size=100):
    poly = nu.polyfit(ra,dec,3)
    ra_array = nu.linspace(minRa, maxRa, size)
    p = nu.poly1d(poly)
    dec_array = p(ra_array)
    return ra_array, dec_array

