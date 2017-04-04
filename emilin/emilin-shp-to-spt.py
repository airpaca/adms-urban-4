#!/usr/bin/env python3
# coding: utf-8

"""Création de données ADMS-Urban 4 (format SPT) des sources routières."""


import configparser
import sys
import pandas
import shapefile
from shapely.geometry import asShape
from emilin import geom


# Arguments
if len(sys.argv) != 2:
    print("usage: adms-shp-to-spt.py configuration.cfg")
    sys.exit(1)
fncfg = sys.argv[1]

# Lecture de la configuration
cfg = configparser.ConfigParser()
cfg.read(fncfg)

fni = cfg.get('fichiers', 'shape')
fno = cfg.get('fichiers', 'output')

field_nom = cfg.get('champs', 'nom')
field_canyon_width = cfg.get('champs', 'canyon_width_m')
field_canyon_height = cfg.get('champs', 'canyon_height_m')
field_height = cfg.get('champs', 'height_m')
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
idx_nom = fields.index(field_nom)
idx_canyon_width = fields.index(field_canyon_width)
idx_canyon_height = fields.index(field_canyon_height)
idx_height = fields.index(field_height)
idx_emissions = {k: fields.index(v) for (k, v) in field_emissions.items()}

# Boucle sur les enregistrements
for row in shp.iterShapeRecords():

    # Split geometry
    g = asShape(row.shape.__geo_interface__)
    lins = geom.split_geom_adms(g)

    for ilin, lin in enumerate(lins, start=1):
        source_name = 'road_{}_{}'.format(row.record[idx_nom], ilin)

        # `spt` file
        canyon_width = row.record[idx_canyon_width]
        canyon_height = row.record[idx_canyon_height]
        height = row.record[idx_height]
        spts.append(
            (source_name, 'No', 0, 0, 'Temperature', 0, 'Actual', 'Velocity', 0,
             0, 'Road', height, 0, canyon_width, canyon_height, 0, 0, 0, 'No',
             2016, 'London (central)', '(Main)', ''))

        # `gpt` file
        gpts.append((source_name, 'Road'))

        # `vgt` file
        for x, y in lin.coords:
            vgts.append((source_name, x, y))

        # `eits` file
        for poladms, idx in idx_emissions.items():
            # kg/an -> g/km/s
            emi = (row.record[idx]
                   * 1000
                   / (3600 * 24 * 365)
                   / (g.length / 1000))
            eits.append((source_name, poladms, emi, 'g/km/s'))

# Entêtes
header_spt = [
    "Source name", "Use VAR file", "Specific heat capacity (J/kg/K)",
    "Molecular mass (g)", "Temperature or density?",
    "Temperature (Degrees C) / Density (kg/m3)",
    "Actual or NTP?", "Efflux type keyword",
    ("Velocity (m/s) / Volume flux (m3/s) / Momentum flux (m4/s2) / "
     "Mass flux(kg/s)"),
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
    f.write("SPTVersion1\n")
    f.write("Traffic dataset name | EFT v6.0.1 (2 VC)\n")
    f.write(spts.to_csv(sep=',', index=False))
print("Sources saved into {fno}.spt".format(**locals()))

# -- gpt
with open(fno + '.gpt', 'w') as f:
    f.write("GPTVersion1\n")
    f.write(gpts.to_csv(sep=',', index=False))
print("Groups saved into {fno}.gpt".format(**locals()))

# -- eit
with open(fno + '.eit', 'w') as f:
    f.write("EITVersion1\n")
    f.write(eits.to_csv(sep=',', index=False))
print("Emissions saved into {fno}.eit".format(**locals()))

# --vgt
with open(fno + '.vgt', 'w') as f:
    f.write("VGTVersion2\n")
    f.write(vgts.to_csv(sep=',', index=False))
print("Geometries saved into {fno}.vgt".format(**locals()))

print("All done !")
