#!/usr/bin/env python3
# coding: utf-8

"""Tests"""

import random
import unittest

from shapely.geometry import LineString

from emilin import geom


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
    def test_split_geom_1(self):
        l1 = random_points(10)
        l2 = random_points(15)
        lt = l1 + l2
        lts = geom.split_geom(lt, [10])
        self.assertEqual(lts[0], l1)
        self.assertEqual(lts[1], l2)

    def test_split_geom_2(self):
        n1, n2, n3 = 3, 5, 4
        l1 = random_points(n1)
        l2 = random_points(n2)
        l3 = random_points(n3)
        lt = l1 + l2 + l3
        lts = geom.split_geom(lt, [n1, n1 + n2])
        self.assertEqual(lts[0], l1)
        self.assertEqual(lts[1], l2)
        self.assertEqual(lts[2], l3)

    def test_split_geom_3(self):
        n1, n2, n3, n4 = 10, 5, 2, 14
        l1 = random_points(n1)
        l2 = random_points(n2)
        l3 = random_points(n3)
        l4 = random_points(n4)
        lt = l1 + l2 + l3 + l4
        lts = geom.split_geom(lt, [n1, n1 + n2, n1 + n2 + n3])
        self.assertEqual(lts[0], l1)
        self.assertEqual(lts[1], l2)
        self.assertEqual(lts[2], l3)
        self.assertEqual(lts[3], l4)

    def test_split_geom_4(self):
        l1 = random_points(5)
        lts = geom.split_geom(l1, [])
        self.assertEqual(lts[0], l1)

    def test_split_geom_adms_1(self):
        l = LineString(random_points(10))
        outs = geom.split_geom_adms(l)
        self.assertEqual(len(outs), 1)
        self.assertEqualGeom(outs[0], l)

    def test_split_geom_adms_2(self):
        l = LineString(random_points(50))
        outs = geom.split_geom_adms(l)
        self.assertEqual(len(outs), 1)
        self.assertEqualGeom(outs[0], l)

    def test_split_geom_adms_3(self):
        p1 = random_points(50)
        p2 = random_points(10)
        l = LineString(p1 + p2)
        outs = geom.split_geom_adms(l)
        self.assertEqual(len(outs), 2)
        self.assertEqualGeom(outs[0], LineString(p1))
        self.assertEqualGeom(outs[1], LineString(p2))

    def test_split_geom_adms_4(self):
        p1 = random_points(50)
        p2 = random_points(50)
        p3 = random_points(20)
        l = LineString(p1 + p2 + p3)
        outs = geom.split_geom_adms(l)
        self.assertEqual(len(outs), 3)
        self.assertEqualGeom(outs[0], LineString(p1))
        self.assertEqualGeom(outs[1], LineString(p2))
        self.assertEqualGeom(outs[2], LineString(p3))


if __name__ == '__main__':
    unittest.main()
