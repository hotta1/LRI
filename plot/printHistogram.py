import sys
import xml.etree.ElementTree
for file in sys.argv[1:]:

  print "plot '-' w l"

  dom = xml.etree.ElementTree.parse(file)

  for p in dom.findall('.//HISTOGRAM/ENTRY'):
#    if p.attrib['name'] == '':
    print p.get('indexvalue'), p.find('VALUE').text

  print "e"
