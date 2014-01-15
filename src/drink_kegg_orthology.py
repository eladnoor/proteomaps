# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 18:42:59 2014

@author: eladn
"""

import urllib, csv, sys, re

KEGG_STEP = 10
ABBREV = 'ko'
FIELD_NAME = 'BRITE'

def FromKeggAPI(s):
    entry_dict = {}
    curr_field = ""
    field_map = {}

    for line in s.split('\n'):
        field = line[0:12].strip()
        value = line[12:]

        if field[:3] == "///":
            bdentry = re.split('\s\s+', field_map['ENTRY'])[0]
            entry_dict[bdentry] = field_map
            field_map = {}
        else:
            if field != "":
                curr_field = field
            if curr_field in field_map:
                field_map[curr_field] = field_map[curr_field] + "\n" + value
            else:
                field_map[curr_field] = value

    if 'ENTRY' in field_map:
        bdentry = re.split('\s\s+', field_map['ENTRY'])[0]
        entry_dict[bdentry] = field_map
    return entry_dict


csv_reader = csv.reader(urllib.urlopen('http://rest.kegg.jp/list/%s/' % ABBREV),
                        delimiter='\t')
all_dbentries = []
for line in csv_reader:
    dbentry = line[0]
    all_dbentries.append(dbentry)
    if len(all_dbentries) > 2:
        break
    
chunks = [all_dbentries[i:i+KEGG_STEP]
          for i in xrange(0, len(all_dbentries), KEGG_STEP)]

for dbentries in chunks:
    sys.stderr.write('Parsing KEGG entities of type %s: %s - %s\n' %
                     (ABBREV, dbentries[0], dbentries[-1]))
   
    url = 'http://rest.kegg.jp/get/' + '+'.join(dbentries)
    sys.stderr.write('%s\n' % url)
    s = urllib.urlopen(url).read()

    entry_dict = FromKeggAPI(s)
    for dbentry, field_map in sorted(entry_dict.iteritems()):
        brite = field_map[FIELD_NAME]
        if brite is not None:
            print dbentry
            print brite
            print '-' * 50
        else:
            sys.stderr.write("WARNING: %s does not contain a BRITE field" %
                             dbentry)
