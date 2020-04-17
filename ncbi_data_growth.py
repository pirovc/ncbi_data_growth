#!/usr/bin/env python3
import pandas as pd
import urllib.request
import re
import matplotlib.pyplot as plt

month_nr = {"Jan":1,"Feb":2,"Mar":3,"Apr":4,"May":5,"Jun":6,"Jul":7,"Aug":8,"Sep":9,"Oct":10,"Nov":11,"Dec":12}
	
def get_file(file,local=False):
	if local:
		with open(file) as f:
			text = "\n".join(f.readlines())
	else:
		response = urllib.request.urlopen(file)
		data = response.read() 
		text = data.decode('utf-8')
	return text

def parse_genbank_wgs(genbank_wgs_file, local=False):
	text = get_file(genbank_wgs_file, local)

	# GenBank
	# start/end anchor
	#Release      Date     Base Pairs   Entries
	#
	#    3    Dec 1982         680338       606
	#   14    Nov 1983        2274029      2427
	#   20    May 1984        3002088      3665	
	# ...
	#   236   Feb 2020   399376854872   216214215
	#   
	#  The following table lists the number of bases and the number of sequence	

	# ? stop greedy regex
	match = re.search('3[ ]+Dec 1982[ ]+680338[ ]+606.*?The', text, re.DOTALL)
	# -1: left overs (The following...)
	# if line: remove extra lines 
	# .strip(): extra spaces before and after ("    3    Dec 1982")
	# .split(): into empty spaces
	# [:5]: get fields of interest
	genbank_entries = [line.strip().split()[:5] for line in match.group(0).split("\n") if line][:-1]

	genbank = pd.DataFrame(genbank_entries, columns=["Release","Month","Year","Base Pairs","Entries"]).dropna()
	genbank['Month'] = genbank['Month'].apply(lambda x: '{0:0>2}'.format(month_nr[x]))
	cols = ["Year","Base Pairs","Entries"]
	genbank[cols] = genbank[cols].apply(pd.to_numeric, errors='ignore', axis=1)
	genbank["Day"] = '{0:0>2}'.format(1)


	# WGS
	# start/end anchor
	#Release      Date     Base Pairs     Entries
	#
	#  129    Apr 2002      692266338      172768
	#  130    Jun 2002     3267608441      397502
	#  131    Aug 2002     3848375582      427771
	# ...
	#  236    Feb 2020  6968991265752  1206720688
	#
	#  The following table provides the number of bases and the number of sequence

	match = re.search('129[ ]+Apr 2002[ ]+692266338[ ]+172768.*?The', text, re.DOTALL)
	wgs_entries = [line.strip().split()[:5] for line in match.group(0).split("\n") if line][:-1]
	
	wgs = pd.DataFrame(wgs_entries, columns=["Release","Month","Year","Base Pairs","Entries"]).dropna()
	wgs['Month'] = wgs['Month'].apply(lambda x: '{0:0>2}'.format(month_nr[x]))
	cols = ["Year","Base Pairs","Entries"]
	wgs[cols] = wgs[cols].apply(pd.to_numeric, errors='ignore', axis=1)
	wgs["Day"] = '{0:0>2}'.format(1)

	return genbank, wgs



def parse_refseq(refseq_file, local=False):

	text = get_file(refseq_file,local)

	# RefSeq
	# start/end anchor
	#3.9 Growth of RefSeq	
	#--------------------
	#Release	Date		Taxons	Nucleotides	 Amino Acids	 Records
	#1	Jun 30, 2003	2005	4672871949	 263588685	 1061675
	#...
	#99	Mar 2, 2020	99842	1865535232080	 64046042055	 231402293
	#
	#Note: Date refers to the data cut-off date, i.e., the release incorporates 
	match = re.search('1[\t]+Jun 30.*?Note', text, re.DOTALL)
	refseq_entries = [line.strip().split()[:8] for line in match.group(0).split("\n") if line][:-1]
	#Nucleotides=Base Pairs
	#Records=Entries
	refseq = pd.DataFrame(refseq_entries, columns=["Release","Month","Day","Year","Taxons","Base Pairs","Amino Acids","Entries"]).dropna()
	refseq["Day"] = refseq["Day"].apply(lambda x: x.rstrip(','))
	refseq['Month'] = refseq['Month'].apply(lambda x: '{0:0>2}'.format(month_nr[x]))
	cols = ["Day","Year","Taxons","Base Pairs","Amino Acids","Entries"]
	refseq[cols] = refseq[cols].apply(pd.to_numeric, errors='ignore', axis=1)
	refseq["Day"] = refseq["Day"].apply(lambda x: '{0:0>2}'.format(x))
	return refseq

#genbank_wgs_file="https://ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt"
genbank_wgs_file = "gbrel.txt"
#refseq_release_nr = get_file("https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER").rstrip()
#refseq_file="https://ftp.ncbi.nlm.nih.gov/refseq/release/release-notes/RefSeq-release" + refseq_release_nr + ".txt"
refseq_file="RefSeq-release99.txt"

genbank, wgs = parse_genbank_wgs(genbank_wgs_file, local=True)
refseq = parse_refseq(refseq_file, local=True)

print(genbank)
print(wgs)
print(refseq)

fig1, ax1 = plt.subplots()
ax1.plot(refseq["Year"], refseq["Base Pairs"])
ax1.plot(wgs["Year"], wgs["Base Pairs"])
plt.show()