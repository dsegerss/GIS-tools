#!/usr/bin/env python
"""
Extract shape-files from OSM data in a Postgis database
"""

from os import path
import sys
from optparse import OptionParser
import subprocess

#Docstrings for the option parser
usage = "usage: %prog [options] "
version = "%prog 1.0"


shpDefList = [{"table": "osm_landusages", "subset": "green",
               "type": ['garden', 'park', 'meadow', 'village_green',
                       'grass', 'golf_course', 'playground', 'pitch',
                       'recreation_ground']},
              {"table": "osm_landusages", "subset": "forest",
               "type": ['wood', 'forest', 'nature_reserve']},
              {"table": "osm_landusages", "subset": "urban",
               "type": ['retail', 'footway', 'pedestrian', 'commercial']},
              {"table": "osm_landusages", "subset": "residential",
               "type": ['residential']},
              {"table": "osm_landusages", "subset": "special_areas",
               "type": ['cemetery', 'college', 'allotments', 'stadium',
                       'sport_centre', 'theatre', 'university', 'school',
                       'place_of_worship', 'library', 'hospital', 'cinema']},
              {"table": "osm_landusages", "subset": "industrial_areas",
               "type": ['industrial', 'quarry']},
              {"table": "osm_landusages", "subset": "farm_areas",
               "type": ['farmland', 'farm', 'farmyard']},
              {"table": "osm_landusages", "subset": "transport_areas1",
               "type": ['parking', 'fuel', 'railway']},
              {"table": "osm_buildings"},
              {"table": "osm_waterareas"},
              {"table": "osm_waterways"},
              {"table": "osm_motorways"},
              {"table": "osm_roads",
               "subset": "with_traffic",
               "type": ["residential", "motorway", "motorway_link",
                        "trunk", "trunk_link", "primary", "primary_link",
                        "secondary", "secondary_link", "tertiary", \
                            "tertiary_link", "unclassified"]},
              {"table": "osm_roads",
               "subset": "residential", "type": ["residential"]},
              {"table": "osm_roads",
               "subset": "motorway", "type": ["motorway", "motorway_link"]},
              {"table": "osm_roads",
               "subset": "trunk", "type": ["trunk", "trunk_link"]},
              {"table": "osm_roads",
               "subset": "primary", "type": ["primary", "primary_link"]},
              {"table": "osm_roads",
               "subset": "secondary", "type": ["secondary", "secondary_link"]},
              {"table": "osm_roads",
               "subset": "tertiary", "type": ["tertiary", "tertiary_link"]},
              {"table": "osm_roads",
               "subset": "unclassified", "type": ["unclassified"]},
              {"table": "osm_roads",
               "subset": "pedestrian",
               "type": ["pedestrian", "footway", "steps",\
                            "cycleway", "bridleway"]},
              {"table": "osm_railways"},
              {"table": "osm_aeroways"},
              {"table": "osm_admin", "subset": "1",
               "admin_level":1},
              {"table": "osm_admin", "subset": "2",
               "admin_level":2},
              {"table": "osm_admin", "subset": "3",
               "admin_level":3},
              {"table": "osm_admin", "subset": "4",
               "admin_level":4},
              {"table": "osm_admin", "subset": "5",
               "admin_level":5},
              {"table": "osm_admin", "subset": "6",
               "admin_level":6},
              {"table": "osm_admin", "subset": "7",
               "admin_level":7},
              {"table": "osm_admin", "subset": "8",
               "admin_level":8},
              {"table": "osm_admin", "subset": "9",
               "admin_level":9},
              {"table": "osm_places", "subset": "locality",
               "type": ['locality']},
              {"table": "osm_places", "subset": "hamlet",
               "type": ['hamlet']},
              {"table": "osm_places", "subset": "village",
               "type": ['village']},
              {"table": "osm_places", "subset": "suburb",
               "type": ['suburb']},
              {"table": "osm_places", "subset": "town",
               "type": ['town']},
              {"table": "osm_places", "subset": "city",
               "type": ['city']},
              {"table": "osm_places", "subset": "region",
               "type": ['region']},
              {"table": "osm_places", "subset": "state",
               "type": ['state']},
              {"table": "osm_places", "subset": "county",
               "type": ['county']}
]


def main():
    #-----------Setting up and unsing option parser-----------------------
    parser = OptionParser(usage=usage, version=version)

    parser.add_option(
        "--doc",
        action="store_true",
        dest="doc",
        help="Prints more detailed documentation and exit"
    )

    parser.add_option(
        "-d", "--dbname",
        action="store",
        dest="dbname",
        default="osm",
        help="Name of database"
    )

    parser.add_option(
        "--clip_table",
        action="store",
        dest="clip_table",
        default="domain_clip_poly",
        help="Name of table containing domain boundaries"
    )

    parser.add_option(
        "-o", "--outdir",
        action="store",
        dest="outdir",
        default=None,
        help="Output directory"
    )

    parser.add_option(
        "-c",
        "--clip",
        action="store",
        dest="clip",
        help="Attribute value to clip"
    )

    parser.add_option(
        "--format",
        action="store",
        dest="format",
        help="Format of extracted files",
        default="ESRI Shapefile"
    )

    parser.add_option(
        "--encoding",
        action="store",
        dest="encoding",
        help="Encoding of extracted files",
        default="HP Roman8"
    )

    parser.add_option(
        "--srs",
        action="store",
        dest="srs",
        help="Proj4 projection definition to be used for extracted files"
    )

    (options, args) = parser.parse_args()

    #------------Process and validate options---------------
    if options.doc:
        print(__doc__)
        return 0

    if len(args) > 0:
        parser.error("Incorrect number of arguments")

    if options.outdir is None:
        parser.error("Need to specify output directory using -o <outdir>")

    for shpDef in shpDefList:
        shpDef["format"] = "\"" + options.format + "\""
        shpDef["clip_table"] = options.clip_table
        shpDef["dbname"] = options.dbname
        shpDef["clip"] = options.clip
        shpDef["srs"] = options.srs
        if "subset" in shpDef:
            shpDef["dest"] = path.join(
                options.outdir, shpDef['table'] + "-" +
                shpDef['subset'] + ".shp")
        else:
            shpDef["dest"] = path.join(
                options.outdir, shpDef['table'] + ".shp")
        query = ("\"SELECT a.* as point FROM " +
                 "%(table)s a, %(clip_table)s b WHERE " % shpDef)
        if "type" in shpDef:
            types = "'" + "','".join(shpDef["type"]) + "'"
            query += " a.type IN (%s) AND " % types
        if "admin_lev" in shpDef:
            query += " a.admin_level= %(admin_level)s AND " % shpDef

        query += " ST_WITHIN(a.geometry,b.geom) "
        query += " AND b.name = '%(clip)s'\"" % shpDef
        shpDef['query'] = query
        shpDef['encoding'] = "\"ENCODING=%s\"" % options.encoding

        if shpDef["srs"] is not None:
            cmd = ("ogr2ogr -f %(format)s " % shpDef +
                   "%(dest)s " % shpDef +
                   "-s_srs \"epsg:3857\" " % shpDef +
                   "-lco %(encoding)s " % shpDef +
                   "-t_srs \"%(srs)s\" " % shpDef +
                   "PG:dbname=%(dbname)s " % shpDef +
                   "-sql %(query)s" % shpDef)
        else:
            cmd = ("ogr2ogr -f %(format)s " % shpDef +
                   "%(dest)s " % shpDef +
                   "-lco %(encoding)s " % shpDef +
                   "PG:dbname=%(dbname)s " % shpDef +
                   "-sql %(query)s" % shpDef)

        if path.exists(shpDef["dest"]):
            print "Extract already exist, remove first"
        else:
            proc = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            retCode = proc.wait()
            if retCode != 0:
                print("Error while running cmd:\n%s" % cmd)
                return 1
            else:
                print("Wrote %(dest)s" % shpDef)

if __name__ == "__main__":
    sys.exit(main())
