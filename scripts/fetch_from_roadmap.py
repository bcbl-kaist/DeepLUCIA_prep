#!/usr/bin/env python

import os
import sys

from urllib import request
from urllib.error import HTTPError
from pathlib import Path

def main():
	if len(sys.argv) < 4:
		print ("Usage : " + sys.argv[0] + " [EID] [marker] [signal_bigwig_filename]")
		sys.exit()
	else:
		EID = sys.argv[1]
		marker = sys.argv[2]
		signal_bigwig_filename = sys.argv[3]
		fetch_from_roadmap(EID,marker,signal_bigwig_filename)

def fetch_from_roadmap(EID,marker,signal_bigwig_filename):
	consolidated_url = "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/{}-{}.pval.signal.bigwig".format(EID,marker)
	imputed_url = "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/{}/{}-{}.imputed.pval.signal.bigwig".format(marker,EID,marker)
	#consolidated_url = f"https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/{EID}-{marker}.pval.signal.bigwig"
	#imputed_url = f"https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/{marker}/{EID}-{marker}.imputed.pval.signal.bigwig"
	Path.mkdir(Path(signal_bigwig_filename).parent, parents=True, exist_ok=True)

	try:
		request.urlopen(consolidated_url)
		request.urlretrieve(consolidated_url,signal_bigwig_filename)
	except HTTPError as error:
		print("consolidated signal is not available. imputed signal will be used instead")
		request.urlretrieve(imputed_url,signal_bigwig_filename)

if __name__ == "__main__":
	main()

