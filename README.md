# Quelques outils pour ADMS-Urban 4


## Création des émissions

Création de fichiers `.spt`, `.gpt`, `.eit`, `.vgt` contenant des informations
sur les sources et leurs émissions.

Ces fichiers peuvent ensuite être importés dans ADMS-Urban 4.

Cf. documentation du modèle ADMS-Urban 4. 


### Sources routières

Script `emilin-shp-to-spt.py`.

Exemple d'un fichier de configuration: `emilin_exemple.cfg`.
__Les émissions doivent être en kg.__

Utilisation:

    python3 emilin-shp-to-spt.py /path/of/configuration.cfg


### Sources cadastrées

Script `cad-shp-to-spt.py`.

Exemple d'un fichier de configuration: `cadastre_exemple.cfg`.
__Les émissions doivent être en kg.__

Utilisation:

    python3 cad-shp-to-spt.py /path/of/configuration.cfg



## Création d'un maillage autour des sources

Script `pymailleur.py`.

Création d'un fichier `.asp` à partir d'un ensemble de fichiers SIG (shapefile) :
 - sources linéaires
 - sources ponctuelles
 - sources surfaciques ou volumiques

Exemple d'un fichier de configuration: `grid_exemple.cfg`.

Utilisation:

    python3 pymailleur.py /path/of/configuration.cfg

