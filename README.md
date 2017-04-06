# Quelques outils pour ADMS-Urban 4


## Création des émissions

Création de fichiers `.spt`, `.gpt`, `.eit`, `.vgt` contenant des informations
sur les sources et leurs émissions.

Ces fichiers peuvent ensuite être importés dans ADMS-Urban 4.

Cf. documentation du modèle ADMS-Urban 4. 


### Sources routières

Script `emilin/emilin_shp_to_spt.py`.

Exemple d'un fichier de configuration: `emilin/emilin_exemple.cfg`.  
__Les émissions doivent être en kg.__

Utilisation:

    python3 emilin/emilin_shp_to_spt.py /path/of/configuration.cfg


### Sources cadastrées

Script `emicad/cad_shp_to_spt.py`.

Exemple d'un fichier de configuration: `emicad/cadastre_exemple.cfg`.  
__Les émissions doivent être en kg.__

Utilisation:

    python3 cad_shp_to_spt.py /path/of/configuration.cfg



## Création d'un maillage autour des sources

Script `mailleur/pymailleur.py`.

Création d'un fichier `.asp` à partir d'un ensemble de fichiers SIG (shapefile) :
 - sources linéaires
 - sources ponctuelles
 - sources surfaciques ou volumiques

Exemple d'un fichier de configuration: `mailleur/grid_exemple.cfg`.

Utilisation:

    python3 mailleur/pymailleur.py /path/of/configuration.cfg

