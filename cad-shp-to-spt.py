#!/usr/bin/env python3
# coding: utf-8

"""Création de données ADMS-Urban 4 (format SPT) du cadastre."""


import sys
import configparser
import shapefile
from shapely.geometry import Polygon
import pandas


# Fonctions
def calc_area(shpobj):
    """Calcul de la surface d'un polygone simple.

    :param shpobj: shapefile._Shape object.
    :return: float.
    """
    assert shpobj.shapeType == shapefile.POLYGON
    return Polygon(shpobj.points).area


# Arguments
if len(sys.argv) != 2:
    print("usage: adms-shp-to-spt-cadastre.py configuration.cfg")
    sys.exit(1)
fncfg = sys.argv[1]

# Lecture de la configuration
cfg = configparser.ConfigParser()
cfg.read(fncfg)

fni = cfg.get('fichiers', 'shape')
fno = cfg.get('fichiers', 'output')

field_id = cfg.get('champs', 'id')
source_height = cfg.getfloat('param', 'hauteur')
txno2 = cfg.getfloat('param', 'txno2')

field_emissions = dict([e.replace(' ', '').split(':')
                        for e in list(cfg['emissions'].values())])

if fno.endswith('.spt'):
    fno = fno[:-4]  # remove ext


# Initialisation des tableaux
spts = []
gpts = []
vgts = []
eits = []


# Ouverture du shapefile
print("Reading data...")
shp = shapefile.Reader(fni)

# Sélection des champs intéressants
fields = [e[0] for e in shp.fields if e[0] != 'DeletionFlag']
idx_nom = fields.index(field_id)
idx_emissions = {k: fields.index(v) for (k, v) in field_emissions.items()}

# Boucle sur les enregistrements
for row in shp.iterShapeRecords():
    source_name = 'grid_' + str(row.record[idx_nom])
    area = calc_area(row.shape)  # compute area of shape in m2

    # `spt` file
    spts.append((source_name, 'No', 0, 0, 'Temperature', 0, 'Actual',
                 'Velocity', 0, 0, 'Grid', 0, 0, source_height, 0, 0, 0, 0,
                 'n/a', 'n/a', 'n/a', '(None)', ''))

    # `gpt` file
    gpts.append((source_name, 'Grid'))

    # `vgt` file
    for x, y in row.shape.points[:4]:  # only the 4 points
        vgts.append((source_name, x, y))

    # `eit` file
    emi_nox = None
    for poladms, idx in idx_emissions.items():

        emi_kg = row.record[idx]
        if emi_kg is None:
            emi_kg = 0

        # kg/an -> g/m2/s
        emi = emi_kg * 1000 / (3600 * 24 * 365) / area
        eits.append((source_name, poladms, emi, 'g/m2/s'))
        if poladms == 'NOx':
            emi_nox = emi
    if emi_nox:
        eits.append((source_name, 'NO2', emi_nox * txno2, 'g/m3/s'))

# Entêtes
header_spt = [
    "Source name", "Use VAR file", "Specific heat capacity (J/kg/K)",
    "Molecular mass (g)", "Temperature or density?",
    "Temperature (Degrees C) / Density (kg/m3)",
    "Actual or NTP?", "Efflux type keyword",
    ("Velocity (m/s) / Volume flux (m3/s) / Momentum flux (m4/s2)"
     " / Mass flux(kg/s)"),
    "Heat release rate (MW)", "Source Type", "Height (m)", "Diameter (m)",
    "Line width (m) / Road width (m) / Volume depth (m) / Grid depth (m)",
    "Canyon height (m)", "Angle 1 (deg)", "Angle 2 (deg)",
    "Mixing ratio (kg/kg)", "Traffic flows used", "Traffic flow year",
    "Traffic flow road type", "Main building", "Comments"]

header_eit = ["Source name", "Pollutant name", "Emission Rate", "Comments"]

header_gpt = ["Source name", "Source Type"]

header_vgt = ["Source name", "X (m)", "Y (m)"]

# Création de dataframe
spts = pandas.DataFrame(spts, columns=header_spt)
gpts = pandas.DataFrame(gpts, columns=header_gpt)
vgts = pandas.DataFrame(vgts, columns=header_vgt)
eits = pandas.DataFrame(eits, columns=header_eit)

# Enregistrement
# -- spt
with open(fno + '.spt', 'w') as f:
    f.write("SPTVersion1\r\n")
    f.write("Traffic dataset name | EFT v6.0.1 (2 VC)\r\n")
    f.write(spts.to_csv(sep=',', index=False, line_terminator='\r\n'))
print("Sources saved into {fno}.spt".format(**locals()))

# -- gpt
with open(fno + '.gpt', 'w') as f:
    f.write("GPTVersion1\r\n")
    f.write(gpts.to_csv(sep=',', index=False, line_terminator='\r\n'))
print("Groups saved into {fno}.gpt".format(**locals()))

# -- eit
with open(fno + '.eit', 'w') as f:
    f.write("EITVersion1\r\n")
    f.write(eits.to_csv(sep=',', index=False, line_terminator='\r\n'))
print("Emissions saved into {fno}.eit".format(**locals()))

# --vgt
with open(fno + '.vgt', 'w') as f:
    f.write("VGTVersion2\r\n")
    f.write(vgts.to_csv(sep=',', index=False, line_terminator='\r\n'))
print("Geometries saved into {fno}.vgt".format(**locals()))

print("All done !")
