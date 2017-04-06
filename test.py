#!/usr/bin/env python3
# coding: utf-8

"""Tests"""

import random
import unittest

from shapely.geometry import LineString

from emilin.emilin_shp_to_spt import split_iterable, split_geom_adms


def random_coord():
    """Random value between 0 and 50."""
    return random.randint(0, 50)


def random_points(n):
    """Create a list of n random points."""
    return [(random_coord(), random_coord()) for i in range(n)]


class ShapelyTestCase(unittest.TestCase):
    """Add assertion function to compare geometries."""

    def assertEqualGeom(self, g1, g2):
        self.assertEqual(g1.type, g2.type)
        self.assertEqual(g1.coords[:], g2.coords[:])


class TestFunctions(ShapelyTestCase):
    def test_split_iterable_1(self):
        iterable = 'abcdefghi'
        idx = []
        objs = split_iterable(iterable, idx)
        self.assertEqual(objs, [iterable])

    def test_split_iterable_2(self):
        iterable = 'abcdefghi'
        idx = [5]
        objs = split_iterable(iterable, idx)
        self.assertEqual(len(objs), 2)
        self.assertEqual(objs[0], 'abcdef')
        self.assertEqual(objs[1], 'fghi')

    def test_split_iterable_3(self):
        iterable = 'abcdefghi'
        idx = [2, 5]
        objs = split_iterable(iterable, idx)
        self.assertEqual(len(objs), 3)
        self.assertEqual(objs[0], 'abc')
        self.assertEqual(objs[1], 'cdef')
        self.assertEqual(objs[2], 'fghi')

    def test_split_iterable_4(self):
        iterable = 'abcdefghi'
        idx = [0, 5]
        objs = split_iterable(iterable, idx)
        self.assertEqual(len(objs), 2)
        self.assertEqual(objs[0], 'abcdef')
        self.assertEqual(objs[1], 'fghi')

    def test_split_iterable_5(self):
        iterable = 'abcdefghi'
        idx = [0, 8]
        objs = split_iterable(iterable, idx)
        self.assertEqual(len(objs), 1)
        self.assertEqual(objs[0], iterable)

    def test_split_iterable_6(self):
        iterable = 'abcdefghi'
        idx = [3, 8]
        objs = split_iterable(iterable, idx)
        self.assertEqual(len(objs), 2)
        self.assertEqual(objs[0], 'abcd')
        self.assertEqual(objs[1], 'defghi')

    def test_split_iterable_points(self):
        l = random_points(20)
        idx = [5, 15]
        lts = split_iterable(l, idx)
        self.assertEqual(len(lts), len(idx) + 1)
        self.assertEqual(lts[0], l[0:6])
        self.assertEqual(lts[1], l[5:16])
        self.assertEqual(lts[2], l[15:])

    def test_split_geom_adms_1(self):
        l = LineString(random_points(10))
        outs = split_geom_adms(l)
        self.assertEqual(len(outs), 1)
        self.assertEqualGeom(outs[0], l)
        for out in outs:
            self.assertLessEqual(len(out.coords[:]), 50)

    def test_split_geom_adms_2(self):
        l = LineString(random_points(50))
        outs = split_geom_adms(l)
        self.assertEqual(len(outs), 1)
        self.assertEqualGeom(outs[0], l)
        for out in outs:
            self.assertLessEqual(len(out.coords[:]), 50)

    def test_split_geom_adms_3(self):
        pts = random_points(60)
        l = LineString(pts)
        outs = split_geom_adms(l)
        self.assertEqual(len(outs), 2)
        self.assertEqualGeom(outs[0], LineString(pts[0:50]))
        self.assertEqualGeom(outs[1], LineString(pts[49:]))
        for out in outs:
            self.assertLessEqual(len(out.coords[:]), 50)

    def test_split_geom_adms_4(self):
        pts = random_points(120)
        l = LineString(pts)
        outs = split_geom_adms(l)
        self.assertEqual(len(outs), 3)
        self.assertEqualGeom(outs[0], LineString(pts[0:50]))
        self.assertEqualGeom(outs[1], LineString(pts[49:99]))
        self.assertEqualGeom(outs[2], LineString(pts[98:]))
        for out in outs:
            self.assertLessEqual(len(out.coords[:]), 50)


if __name__ == '__main__':
    unittest.main()
