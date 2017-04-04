#!/usr/bin/env python3
# coding: utf-8

"""Functions."""


from typing import Union, List
from shapely.geometry import LineString, MultiLineString


def split_geom(points: List, steps: List):
    """Split geometry at spectific points.

    :param points: list of points (list of (x, y) values).
    :param steps: list of index where to cut.
    :return: list of objects.
    """
    if not steps:
        return [points]

    # compute size of each objets from indexes
    lenobjs = [steps[0]] + [e2 - e1 for (e1, e2) in zip(steps, steps[1:])]

    parts = list()
    residue = points[:]
    for lenobj in lenobjs:
        parts.append(residue[0:lenobj])
        residue = residue[lenobj:]
    parts.append(residue)
    return parts


def split_geom_adms(
        geom: Union[LineString, MultiLineString]
        ) -> List[LineString]:
    """Split geometry for ADMS-Urban.
    Separate each multipart and split parts than more of 50 vertex.
    
    :param geom: shapely LineString or MultiLineString.
    :return: list of shapely LineString.
    """
    if geom.type == 'MultiLineString':
        lins = list(geom)
    elif geom.type == 'LineString':
        lins = [geom]
    else:
        raise TypeError("geom type is not LineString or MultiLineString!")

    outs = list()
    for lin in lins:
        steps = list(range(50, len(lin.coords), 50))
        splits = split_geom(lin.coords[:], steps)
        for split in splits:
            outs.append(LineString(split))
    return outs
