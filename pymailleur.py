#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-

""" Création d'un maillage personnalisé pour les modèles urbains. """


import os
import configparser
import logging
from math import cos, sin
import numpy
import shapefile
import sys
from shapely.geometry import Point, LineString, asShape
from shapely.ops import unary_union


class Pts:
    """ Point with name and source type.
    """
    def __init__(self, name, typsrc, x, y):
        self.name = name
        self.typsrc = typsrc
        self.p = Point(x, y)

    @property
    def x(self):
        return self.p.x

    @property
    def y(self):
        return self.p.y

    @property
    def xy(self):
        return self.p.x, self.p.y

    @property
    def wkt(self):
        return self.p.wkt

    def __str__(self):
        return "Pts('%s', %s, (%s, %s))" % (self.name, self.typsrc, self.x,
                                            self.y)

    def __repr__(self):
        return "<%s>" % self.__str__()


def shift(p, dx=0, dy=0):
    """ Move point.
    :param p: shapely Point object or Pts object.
    :param dx: X distance.
    :param dy: Y distance.
    :return: shapely Point object.
    """
    if isinstance(p, Point):
        return Point(p.x + dx, p.y + dy)
    elif isinstance(p, Pts):
        return Pts(p.name, p.typsrc, p.x + dx, p.y + dy)
    else:
        raise SystemError("cannot reconize 'p' object")


def remove_duplicate_points(l):
    """ Remove duplicate points (function of x and y only)
    :param l: list of Pts object.
    :return: list of Pts object.
    """
    l2, xy = [], []
    for p in l:
        if p.xy not in xy:
            xy.append(p.xy)
            l2.append(p)
    del xy
    return l2


class Mailleur:
    """ Création d'un maillage personnalisé pour les modèles urbains. """

    def __init__(self, xmin, xmax, ymin, ymax, logname='main', pfx_log=""):
        """ Mailleur sur une zone.
        """
        if xmin > xmax:
            raise ValueError("xmin > xmax")
        if ymin > ymax:
            raise ValueError("ymin > ymax")
        self.log = logging.getLogger('%s.pymailleur' % logname)
        self.pfx_log = pfx_log
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.ptsreg = []  # grille régulière
        self.ptsponc = []  # grille autour des points
        self.ptslin = []  # grille autour des linéaires
        self.ptssurf = []  # grille autour des surfaciques

    def __len__(self):
        """ Nombre de points. """
        return (len(self.ptslin) + len(self.ptsponc) + len(self.ptsreg) +
                len(self.ptssurf))

    def edge(self, loc):
        """ Retourne la limite d'un côté du domaine.
        :param loc: char représentant le nom du côté
                    - l:left, r:right, t:top, b:bottom
        :return: LineString shapely object.
        """
        if loc == 'l':  # left
            p1, p2 = (self.xmin, self.ymin), (self.xmin, self.ymax)
        elif loc == 'r':  # right
            p1, p2 = (self.xmax, self.ymin), (self.xmax, self.ymax)
        elif loc == 't':  # top
            p1, p2 = (self.xmin, self.ymax), (self.xmax, self.ymax)
        elif loc == 'b':  # bottom
            p1, p2 = (self.xmin, self.ymin), (self.xmax, self.ymin)
        else:
            return None
        return LineString([p1, p2])

    @staticmethod
    def __splitmulti(geom):
        """ Séparation de Multi* en *.
        :param geom: a shapely object.
        """
        simplegeoms = ['Point', 'LineString', 'LinearRing', 'Polygon']
        multigeoms = ['MultiPoint', 'MultiLineString', 'MultiPolygon']
        if geom.geom_type in simplegeoms:
            return [geom, ]
        elif geom.geom_type in multigeoms:
            return list(geom)
        else:
            raise ValueError(geom.geom_type)

    @staticmethod
    def __points_on_lines(lins, dx):
        """ Création de points régulièrements espacés sur les lignes.
        :param lins: liste de lignes.
        :param dx: espacement entre les points.
        """
        pts = []
        for lin in lins:
            for i in range(0, int(lin.length / dx) * dx + dx, dx):
                # Interpolation et coords du point
                [(x, y)] = lin.interpolate(i).coords[:]
                pts.append(Point(x, y))
        return pts

    def add_regular(self, dx, dy, dlx=None, dly=None):
        """ Création d'un maillage régulier sur la zone.
            Ajout de points le long des limites du domaine.
        :param dx: résolution X de la grille.
        :param dy: résolution Y de la grille.
        :param dlx: résolution X sur la limite du domaine.
        :param dly: résolution Y sur la limite du domaine.
        """

        if dlx is None:
            dlx = dx
        if dly is None:
            dly = dy

        # -- 1. grille régulière
        lx = numpy.arange(self.xmin + dx, self.xmax, dx)
        ly = numpy.arange(self.ymin + dy, self.ymax, dy)
        i = 1
        for x in lx:
            for y in ly:
                pts = Pts('reg_c_%i' % i, 'reg', x, y)
                self.ptsreg.append(pts)
                del pts
                i += 1
        del lx, ly

        # -- 2. points supplémentaires le long des limites du domaine
        i = 1
        for x in [self.xmin, self.xmax]:
            for y in numpy.arange(self.ymin, self.ymax + dly, dly):
                pts = Pts('reg_lim_%i' % i, 'reg', x, y)
                self.ptsreg.append(pts)
                del pts
                i += 1
        for y in [self.ymin, self.ymax]:
            for x in numpy.arange(self.xmin, self.xmax + dlx, dlx):
                pts = Pts('reg_lim_%i' % i, 'reg', x, y)
                self.ptsreg.append(pts)
                del pts
                i += 1

        # -- 3. remove duplicate
        self.ptsreg = remove_duplicate_points(self.ptsreg)
        self.log.info("{pfx} Création de {n} points réguliers".format(
            pfx=self.pfx_log, n=len(self.ptsreg)
        ))

    def add_ponct(self, points, r):
        """ Création d'un maillage autour d'une liste de points.
        :param points: list of shapely.geometry.Point object.
        :param r: rayon.
        """
        dist = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]  # % du rayon
        for ipt, p in enumerate(points, start=1):

            # Add source point
            pts = Pts('ponct_p%i' % ipt, 'ponct', p.x, p.y)
            self.ptsponc.append(pts)

            for ith, theta in enumerate(numpy.arange(0, 360, 18), start=1):

                # Coords des points sur le cercle de rayon r
                xc, yc = p.x + r * cos(theta), p.y + r * sin(theta)

                # Ligne entre le centre (source) et le cercle
                l = LineString([(p.x, p.y), (xc, yc)])

                for idist, pc in enumerate(dist, start=1):

                    # Interpolation et coords du point
                    [(x, y)] = l.interpolate(pc, normalized=True).coords[:]
                    pts = Pts('ponct_p%i_th%i_d%i' % (ipt, ith, idist),
                              'ponct', x, y)
                    self.ptsponc.append(pts)

        # Remove duplicate
        self.ptsponc = remove_duplicate_points(self.ptsponc)

        self.log.info(
            "{pfx} Création de {n} points autour des sources "
            "ponctuelles".format(
                pfx=self.pfx_log, n=len(self.ptsponc)))

    def add_lin_zt(self, lins, param):
        """ Création d'un maillage autour de linéaire avec la méthode des
            zones tampons.
            Ajout de points aux extrémités du domaine.
        :param lins: list of shapely.geometry.LineString or MultiLineString.
        :param param: list of (dzt, dx).
        """
        # 24/07/2012 -- bug (?) on ulin.buffer() --> create null geometry
        # union(buffer of each linestring) instead of
        # buffer(union of linestring)
        for (dzt, dx) in param:

            # -- 1. points sur et autour des axes
            if dzt == 0:
                for ilin, lin in enumerate(lins, start=1):
                    boundaries = self.__splitmulti(lin)
                    pts_on_lines = self.__points_on_lines(boundaries, dx)
                    for ipt, pt in enumerate(pts_on_lines, start=1):
                        pts = Pts('lin_dzt%s_lin%i_%i' % (dzt, ilin, ipt),
                                  'lin', pt.x, pt.y)
                        self.ptslin.append(pts)
                        del pts
                    del boundaries, pts_on_lines
            else:
                # Buffer for each linestring
                bufs = [l.buffer(dzt) for l in lins]
                uzts = self.__splitmulti(unary_union(bufs))  # union of buffers
                for iuzt, uzt in enumerate(uzts, start=1):
                    boundaries = self.__splitmulti(uzt.boundary)
                    pts_on_lines = self.__points_on_lines(boundaries, dx)
                    for ipt, pt in enumerate(pts_on_lines, start=1):
                        pts = Pts('lin_dzt%s_uzt%i_%i' % (dzt, iuzt, ipt),
                                  'lin', pt.x, pt.y)
                        self.ptslin.append(pts)
                        del pts
                    del boundaries, pts_on_lines
                del uzts, bufs
            self.log.debug(
                "add_lin_zt: points, dzt = {dzt}, dx = {dx}".format(
                    dzt=dzt, dx=dx))

            # -- 2. points aux limites du domaine
            for loc in ['l', 'r', 't', 'b']:
                edge = self.edge(loc)
                i = 1

                # Intersection linéaire / limites du domaine
                ptints = [edge.intersection(lin)
                          for lin in lins if edge.intersects(lin)]
                for pis in ptints:
                    for pi in self.__splitmulti(pis):  # split multipoints
                        self.ptslin.append(Pts('lin_edge-%s_%i' % (loc, i),
                                               'lin', pi.x, pi.y))
                        i += 1
                        if loc in ['l', 'r']:
                            for psh in [shift(pi, dy=dzt), shift(pi, dy=-dzt)]:
                                self.ptslin.append(
                                    Pts('lin_edge-%s_%i' % (loc, i),
                                        'lin', psh.x, psh.y))
                                i += 1
                        if loc in ['t', 'b']:
                            for psh in [shift(pi, dx=dzt), shift(pi, dx=-dzt)]:
                                self.ptslin.append(
                                    Pts('lin_edge-%s_%i' % (loc, i),
                                        'lin', psh.x, psh.y))
                                i += 1
            self.log.debug(
                "add_lin_zt: points aux limites, dzt = {dzt}".format(dzt=dzt))

        del lins  # , ulin
        self.log.info(
            "{pfx} Création de {n} points autour des sources linéaires".format(
                pfx=self.pfx_log, n=len(self.ptslin)))

    def add_surf(self, surfs, param_reg, param_zt):
        """ Création d'un maillage régulier dans le polygone et d'un maillage
            autour de polygones avec la méthode des zones tampons.
        :param surfs: list of shapely.geometry.Polygon or MultiPolygon object.
        :param param_reg: regular grid resolution.
        :param param_zt: list of (dzt, dx).
        """
        # Fusion et désagrégation des polygones
        for isurf, surf in enumerate(self.__splitmulti(unary_union(surfs)),
                                     start=1):
            for poly in self.__splitmulti(surf):  # split multipolygon

                # -- 1. grille régulière dans le polygone
                xmin, ymin, xmax, ymax = poly.bounds  # emprise
                lx = numpy.arange(xmin, xmax + param_reg, param_reg)
                ly = numpy.arange(ymin, ymax + param_reg, param_reg)
                i = 1
                for x in lx:
                    for y in ly:
                        if poly.contains(Point(x, y)):
                            self.ptssurf.append(
                                Pts('surf_%i_in%i' % (isurf, i), 'surf', x, y))
                            i += 1

                # -- 2. grille autour du polygone
                for (dzt, dx) in param_zt:
                    if dzt == 0:
                        exterior = LineString(poly.exterior.coords)
                        pts_surf = self.__points_on_lines([exterior, ], dx)
                        for iptsrf, ptsrf in enumerate(pts_surf, start=1):
                            pts = Pts('surf_%i_dzt%f_%i' % (isurf, dzt, iptsrf),
                                      'surf', ptsrf.x, ptsrf.y)
                            self.ptssurf.append(pts)
                            del pts
                    else:
                        zt = poly.buffer(dzt)
                        exterior = LineString(zt.exterior.coords)
                        pts_surf = self.__points_on_lines([exterior, ], dx)
                        for iptsrf, ptsrf in enumerate(pts_surf, start=1):
                            pts = Pts('surf_%i_dzt%f_%i' % (isurf, dzt, iptsrf),
                                      'surf', ptsrf.x, ptsrf.y)
                            self.ptssurf.append(pts)
                            del pts

        self.log.info(
            "{pfx} Création de {n} points pour les sources "
            "surfaciques / volumiques".format(
                pfx=self.pfx_log, n=len(self.ptssurf)))

    def export_asp(self, fn, z=1.5, fmt=".0f"):
        """ Exportation des points au format ASP pour ADMS-Urban.
        :param fn: filename.
        :param z: hauteur Z des points.
        :param fmt: format des coordonnées X et Y.
        """
        f = open(fn, 'w')
        for pts in self.ptsreg + self.ptsponc + self.ptslin + self.ptssurf:
            if self.xmin <= pts.x <= self.xmax and \
                                    self.ymin <= pts.y <= self.ymax:
                f.write("{nom},{x},{y},{z}\n".format(
                    nom=pts.name,
                    x=format(pts.x, fmt),
                    y=format(pts.y, fmt),
                    z=z))
        f.close()
        self.log.info("{pfx} Export asp de {n} points ok".format(
            pfx=self.pfx_log, n=len(self)))

    def export_shp(self, fn, z=1.5):
        """ Exportation des points au format Shapefile.
        :param fn: filename.
        :param z: hauteur Z des points.
        """
        w = shapefile.Writer(shapeType=shapefile.POINT)
        w.field('oid', 'N', 11, 0)
        w.field('name', 'C')
        w.field('typsrc', 'C')
        w.field('x', 'N', 18, 1)
        w.field('y', 'N', 18, 1)
        w.field('z', 'N', 18, 1)
        i = 0
        for pts in self.ptsreg + self.ptsponc + self.ptslin + self.ptssurf:
            if self.xmin <= pts.x <= self.xmax and \
                                    self.ymin <= pts.y <= self.ymax:
                w.record(i, pts.name, pts.typsrc, pts.x, pts.y, z)
                w.point(pts.x, pts.y)
                i += 1
        w.save(fn)
        self.log.info("{pfx} Export shapefile de {n} points ok".format(
            pfx=self.pfx_log, n=len(self)))


if __name__ == '__main__':

    # Log
    log = logging.getLogger('pymailleur')
    log.setLevel(logging.INFO)

    lch = logging.StreamHandler()
    log.addHandler(lch)

    # Intro
    log.info("Création d'un maillage...")
    log.info("Lecture de la configuration")

    # Lecture de la configuration
    fncfg = sys.argv[1]
    cfg = configparser.ConfigParser()
    cfg.read(fncfg)

    if 'output' in cfg:
        dirout = cfg.get('output', 'dir')
    else:
        dirout = '.'

    cfg_reg_shp = cfg.get('regular', 'shp')
    cfg_reg_dx = cfg.getint('regular', 'dx')
    cfg_reg_dy = cfg.getint('regular', 'dy')
    log.info(" - maillage régulier: {dx} m x {dy} m".format(
        dx=cfg_reg_dx, dy=cfg_reg_dy))

    cfg_lin_shp = cfg.get('srclin', 'shp')
    cfg_lin_zt = cfg.get('srclin', 'zt')
    if cfg_lin_shp and cfg_lin_zt:
        cfg_lin_zt = tuple(eval(cfg_lin_zt))
        log.info((" - maillage irrégulier autour des sources linéaires: "
                  "{cfg_lin_zt} (m, m)").format(**locals()))

    cfg_ponct_shp = cfg.get('srcponct', 'shp')
    cfg_ponct_rayon = cfg.get('srcponct', 'rayon')
    if cfg_ponct_shp and cfg_ponct_rayon:
        cfg_ponct_rayon = int(cfg_ponct_rayon)
        log.info((" - maillage irrégulier autour des sources ponctuelles: "
                  "rayon de {cfg_ponct_rayon} m").format(**locals()))

    cfg_surf_shp = cfg.get('srcsurf', 'shp')
    cfg_surf_reg = cfg.get('srcsurf', 'reg')
    cfg_surf_zt = cfg.get('srcsurf', 'zt')
    if cfg_surf_shp and cfg_surf_reg and cfg_surf_zt:
        cfg_surf_reg = int(cfg_surf_reg)
        cfg_surf_zt = tuple(eval(cfg_surf_zt))
        log.info((" - maillage irrégulier autour des sources surfaciques: "
                  "reg {cfg_surf_reg} m, {cfg_surf_zt} (m, m)").format(
            **locals()))

    # Emprise du domaine
    shpreg = shapefile.Reader(cfg_reg_shp)
    log.info("Calcul de l'étendue géographique")
    xmin = min([obj.bbox[0] for obj in shpreg.shapes()])
    ymin = min([obj.bbox[1] for obj in shpreg.shapes()])
    xmax = max([obj.bbox[2] for obj in shpreg.shapes()])
    ymax = max([obj.bbox[3] for obj in shpreg.shapes()])
    log.info(" - {xmin:.0f} m <= x <= {xmax:.0f} m  ({lx:.0f} m)".format(
        lx=(xmax - xmin), **locals()))
    log.info(" - {ymin:.0f} m <= y <= {ymax:.0f} m  ({ly:.0f} m)".format(
        ly=(ymax - ymin), **locals()))

    # Open files
    log.info("Lecture de la géométrie des sources...")
    if cfg_lin_shp and cfg_lin_zt:
        shplin = shapefile.Reader(cfg_lin_shp)
        lins = [asShape(s.__geo_interface__) for s in shplin.shapes()]
        log.info(" - {} sources linéaires".format(len(lins)))
    else:
        lins = None

    if cfg_ponct_shp and cfg_ponct_rayon:
        shpponct = shapefile.Reader(cfg_ponct_shp)
        poncts = [asShape(s.__geo_interface__) for s in shpponct.shapes()]
        log.info(" - {} sources ponctuelles".format(len(poncts)))
    else:
        poncts = None

    if cfg_surf_shp and cfg_surf_reg and cfg_surf_zt:
        shpsurf = shapefile.Reader(cfg_surf_shp)
        surfs = [asShape(s.__geo_interface__) for s in shpsurf.shapes()]
        log.info(" - {} sources surfaciques/volumiques".format(len(surfs)))
    else:
        surfs = None

    # Création du maillage
    log.info("Création du maillage...")
    m = Mailleur(xmin, xmax, ymin, ymax, logname='pymailleur', pfx_log=" -")
    m.add_regular(cfg_reg_dx, cfg_reg_dy)
    if lins and cfg_lin_zt:
        m.add_lin_zt(lins, cfg_lin_zt)
    if poncts and cfg_ponct_rayon:
        m.add_ponct(poncts, cfg_ponct_rayon)
    if surfs and cfg_surf_reg and cfg_surf_zt:
        m.add_surf(surfs, cfg_surf_reg, cfg_surf_zt)

    # Export du résultat
    log.info("Export des résultats...")
    bn = os.path.basename(fncfg)[:-4]
    m.export_asp(os.path.join(dirout, '{}.asp'.format(bn)), fmt=".1f")
    log.info("   - {}/{}.asp".format(dirout, bn))
    m.export_shp(os.path.join(dirout, '{}.shp'.format(bn)))
    log.info("   - {}/{}.shp".format(dirout, bn))
