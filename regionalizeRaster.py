#!/usr/bin/env python
# -*- coding: latin-1 -*-
"""
****************************************************************************
Name: regionalizeRaster.py
Created: 30 May 2013
Author: David Segersson

Description
----------------------------------------------------------------------------
Redistribute a raster so that regions are scaled to match a statistics table.
1. The sum for each region is calculated
2. All cells in the region are divided by this sum (normalized)
3. The cells in the region are multiplied with the corresponding sum found
in the statistics table.
"""

# standard modules
import sys
from optparse import OptionParser
import logging

# third party modules
import numpy as np

# pyAirviro
from pyAirviro.geo.raster import Raster
from pyAirviro.other.logging import get_loglevel
from pyAirviro.other.datatable import DataTable

# ------Global variables--------
log = None
# -----------------------------


def distributeRegion(argList):
    regionArray = argList[0]
    data = argList[1]
    code = argList[2]
    result = argList[3]
    regionalTotal = argList[4]

    mask = np.array(regionArray == code)
    regKey = mask * data
    regSum = np.sum(np.sum(regKey))
    # These two are put here to ensure that all threads are executed
    # equally fast.
    onesRast = mask * (np.ones(mask.shape) * 1.0)
    regMaskSum = np.sum(np.sum(onesRast))
    if regSum > 0:
        result = result + regKey / regSum
    elif regionalTotal > 0:
        log.warning(
            "For code: %i the raster is zero but regional total is: %f," % (
                code, regionalTotal) +
            "the distribution is homogenious over the whole region")
        result = result + onesRast / regMaskSum
    return result

# Docstrings for the option parser
usage = "usage: %prog [options] "
version = "%prog 1.0"


def main():
    # -----------Setting up and unsing option parser-----------------------
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("-d", "--doc",
                      action="store_true", dest="doc",
                      help="Prints more detailed documentation and exit")

    parser.add_option("-v",
                      action="store_const", const=logging.DEBUG,
                      dest="loglevel", default=get_loglevel(),
                      help="Produce verbose output")

    parser.add_option("-o", "--output",
                      action="store", dest="outfileName", default=None,
                      help="Output file")

    parser.add_option("-i", "--i",
                      action="store", dest="infile",
                      help="Raster to be regionalized")

    parser.add_option("-r", "--regdef",
                      action="store", dest="regdef",
                      help="Raster defining the regions found " +
                      "in statistics table")

    parser.add_option("-s", "--stat",
                      action="store", dest="stat",
                      help="Tab-separated statistics table, " +
                      "with 1st column representing geocode and value " +
                      "representing region sum")

    parser.add_option("-c", "--colname",
                      action="store", dest="colname",
                      help="Header of column in stat table to fetch " +
                      "region sums from")

    (options, args) = parser.parse_args()

    # ------------Setting up logging capabilities -----------
    rootLogger = logger.RootLogger(int(options.loglevel))
    global log
    log = rootLogger.getLogger(sys.argv[0])

    # ------------Process and validate options---------------
    if options.doc:
        print __doc__
        sys.exit()

    if len(args) > 0:
        parser.error("Incorrect number of arguments")

    regdef = Raster()
    if options.regdef is not None:
        regdef.read(options.regdef)
    else:
        log.error("No region definition raster specified")
        sys.exit(1)

    key = Raster()
    if options.infile is not None:
        key.read(options.infile)
    else:
        log.info("No initial distribution raster given, using " +
                 "homogenious distribution")
        key.assign(regdef)
        key.data = np.where(regdef.data != regdef.nodata, 1, regdef.nodata)

    if options.stat is None:
        log.error("No statistics table specified")
        sys.exit(1)

    if (regdef.xll != key.xll or regdef.yll != key.yll):
        log.error("The raster domain size is differs between key raster " +
                  "and region raster")

    stat = DataTable()
    stat.read(options.stat, defaultType=float)
    print stat.desc
    stat.convertCol(stat.desc[0]["id"], int)
    print stat.desc
    # Set first column as unique id for rows
    stat.setKeys([stat.desc[0]["id"]])

    # Remove nodata in rasters
    key.nodataToZero()
    regdef.nodataToZero()

    log.info("Consistency och completeness check for codes " +
             "in regdef and statistics")
    # Create list of codes in raster
    regdefCodes = regdef.unique()
    regdefCodes = map(int, regdefCodes)

    if 0.0 in regdefCodes:
        regdefCodes.remove(0.0)

    statRegCodes = [row[0] for row in stat.data]

    if options.colname is not None:
        try:
            col = stat.colIndex[options.colname]
        except KeyError:
            log.error("No column named " +
                      "%s found in statistics table" % options.colname)
            sys.exit(1)
    else:
        col = 1
    colId = stat.desc[col]["id"]

    errFound = False
    for code in statRegCodes:
        # Assuring that the regional codes are present in regdef
        if code not in regdefCodes:
            log.error("Input Error: code:" + str(int(code)) +
                     "in reg totals is not represented in regional raster")
            errFound = True

    # Assuring that the regional raster IDs are present
    # in the regional statistics
    for code in regdefCodes:
        if code not in statRegCodes:
            log.error("Input Error: ID:" + str(code) +
                     " in region raster is not represented in " +
                      "regional totals!")
            errFound = True

    if errFound:
        sys.exit(1)

    # For all regions, calculate the regional sum,
    # set the key values to key/regSum
    res = np.zeros(key.data.shape)
    log.info("Regionalizing key raster")
    for code in regdefCodes:
        emis = stat.lookup(colId, code)
        log.debug("Statistics for reg %i: %f" % (code, emis))
        mask = np.array(regdef.data == code)
        regKey = mask * key.data
        regSum = np.sum(np.sum(regKey))
        if regSum > 0:
            res = res + emis * regKey / regSum
        else:
            log.warning(
                "Distribution key is zero for geocode" +
                " %i, using homogenuous distribution for this region" % code)
            res = res + emis * mask / np.sum(np.sum(mask))

    key.data = res

    # resultRast=(totFractionRast*key)*regionalTotals.regSum(substance)
    key.write(options.outfileName)
    log.info("Result with sum: %f written" % res.sum())

if __name__ == "__main__":
    main()
