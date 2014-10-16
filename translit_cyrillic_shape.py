#!/usr/bin/env python
"""
Translit cyrillic fields in UTF-8 shape-file
"""
from os import path
import sys
from optparse import OptionParser

try:
    from osgeo import ogr
except ImportError:
    print("Could not load ogr. Must have GDAL installed")
    sys.exit(1)


#Docstrings for the option parser
usage = "usage: %prog [options] "
version = "%prog 1.0"

cyrillic_translit = {
    u'\u0410': 'A', u'\u0430': 'a',
    u'\u0411': 'B', u'\u0431': 'b',
    u'\u0412': 'V', u'\u0432': 'v',
    u'\u0413': 'G', u'\u0433': 'g',
    u'\u0414': 'D', u'\u0434': 'd',
    u'\u0415': 'E', u'\u0435': 'e',
    u'\u0416': 'Zh', u'\u0436': 'zh',
    u'\u0417': 'Z', u'\u0437': 'z',
    u'\u0418': 'I', u'\u0438': 'i',
    u'\u0419': 'I', u'\u0439': 'i',
    u'\u041a': 'K', u'\u043a': 'k',
    u'\u041b': 'L', u'\u043b': 'l',
    u'\u041c': 'M', u'\u043c': 'm',
    u'\u041d': 'N', u'\u043d': 'n',
    u'\u041e': 'O', u'\u043e': 'o',
    u'\u041f': 'P', u'\u043f': 'p',
    u'\u0420': 'R', u'\u0440': 'r',
    u'\u0421': 'S', u'\u0441': 's',
    u'\u0422': 'T', u'\u0442': 't',
    u'\u0423': 'U', u'\u0443': 'u',
    u'\u0424': 'F', u'\u0444': 'f',
    u'\u0425': 'Kh', u'\u0445': 'kh',
    u'\u0426': 'Ts', u'\u0446': 'ts',
    u'\u0427': 'Ch', u'\u0447': 'ch',
    u'\u0428': 'Sh', u'\u0448': 'sh',
    u'\u0429': 'Shch', u'\u0449': 'shch',
    u'\u042a': '"', u'\u044a': '"',
    u'\u042b': 'Y', u'\u044b': 'y',
    u'\u042c': "'", u'\u044c': "'",
    u'\u042d': 'E', u'\u044d': 'e',
    u'\u042e': 'Iu', u'\u044e': 'iu',
    u'\u042f': 'Ia', u'\u044f': 'ia'
    }


def transliterate(word, translit_table):
    converted_word = ''
    for char in word:
        transchar = ''
        if char in translit_table:
            transchar = translit_table[char]
        else:
            transchar = char
        converted_word += transchar
    return converted_word


def main():
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-i",
                      action="store", dest="infile", default=None,
                      help="Shapefile to update")

    parser.add_option("-f", '--fields',
                      action="store", dest="field_list",
                      help="Comma-separated list of fields to transliterate")

    (options, args) = parser.parse_args()

    if options.infile is None:
        print("No infile specified")
        sys.exit(0)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.Open(options.infile, 1)
    if data_source is None:
        print("Could not open shape file")
        sys.exit(1)
    layer = data_source.GetLayer()
    feature = layer.GetNextFeature()
    if options.field_list is not None:
        field_names = options.field_list.split(',')
    else:
        print("No field names specified")
        sys.exit(0)

    while feature:
        for fname in field_names:
            try:
                val = feature.GetFieldAsString(fname.encode("latin-1"))
            except:
                print("Warning: could not find field %s" % fname)
                field_names.remove(fname)
                continue
            trans_val = transliterate(val.decode("UTF-8"), cyrillic_translit)
            feature.SetField(fname, trans_val.encode("UTF-8"))
        layer.SetFeature(feature)
        feature = layer.GetNextFeature()
    data_source.Destroy()

if __name__ == "__main__":
    main()
