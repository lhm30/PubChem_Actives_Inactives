import cPickle
import pubchempy
from pubchempy import get_compounds
import Bio
from Bio import Entrez
import StringIO
import urllib
import urllib2
import socket
import sys
import os
import time
import re
import gzip
import numpy as np
from optparse import OptionParser
from functools import wraps
from lxml import etree as ET
import re

#optionparser options
parser = OptionParser()
parser.add_option('-i','--input', dest='inf', help='Input file (Uniprots [default] or EGIDs [must supply the --geneids option])', metavar='FILE')
parser.add_option('--geneids', action='store_true', default=False, dest='egids', help='Toggle Entrez Gene ID input format')
parser.add_option('--actives', action='store_true', default=False, dest='actives', help='Toggle extract actives')
parser.add_option('--inactives', action='store_true', default=False, dest='inactives', help='Toggle extract inactives')
(options, args) = parser.parse_args()


def check_email_address(address):
  is_valid = re.search('[^@]+@[^@]+\.[^@]+', address)
  if is_valid: return True
  else:
    print 'Invalid email'
    return False

msg = 'Please enter your e-mail address for the NCBI Eutils: '
email_address = raw_input(msg)
while not check_email_address(email_address):
  # Keep prompting for email if not valid
  email_address = raw_input(msg)
else: pass

Entrez.email = email_address
socket.setdefaulttimeout(15)
opener = urllib2.build_opener()
opener.addheaders = [('User-agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_8_4) AppleWebKit/536.30.1 (KHTML, like Gecko) Version/6.0.5 Safari/536.30.1')]
opener.addheaders = [('Connection', 'keep-alive')]
opener.addheaders = [('Accept', 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8')]
opener.addheaders = [('Accept-Language', 'en-us')]


def introMessage():
	print '=============================================================================================='
	print ' Get PubChem Inactives or Actives'
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return

def retry(ExceptionToCheck, tries=20, delay=4, backoff=2, logger=None):
	def deco_retry(f):
		@wraps(f)
		def f_retry(*args, **kwargs):
			mtries, mdelay = tries, delay
			while mtries > 1:
				try:
					return f(*args, **kwargs)
				except ExceptionToCheck, e:
					msg = "\n" + time.ctime() + " %s, Retrying in %d seconds..." % (str(e), mdelay)
					if logger:
						print msg
						logger.warning(msg)
					else:
						print msg
						file = open('logger.txt', 'a')
						file.write(msg)
						file.close()
					mtries -= 1
					mdelay *= backoff
			return f(*args, **kwargs)
		return f_retry
	return deco_retry

def retry2(ExceptionToCheck, tries=3, delay=2, backoff=2, logger=None):
	def deco_retry(f):
		@wraps(f)
		def f_retry(*args, **kwargs):
			mtries, mdelay = tries, delay
			while mtries > 1:
				try:
					return f(*args, **kwargs)
				except ExceptionToCheck, e:
					msg = "\n" + time.ctime() + " %s, Retrying in %d seconds..." % (str(e), mdelay)
					if logger:
						logger.warning(msg)
					else:
						file = open('logger2.txt', 'a')
						file.write(msg)
						file.close()
					mtries -= 1
					mdelay *= backoff
			return f(*args, **kwargs)
		return f_retry
	return deco_retry
	
@retry(Exception)
def esearch_retry(inp1,inp2):
	record = Entrez.esearch(db=inp1, term=inp2, retmax = 99999)
	return Entrez.read(record)

@retry(Exception)
def elink_retry(inp1,inp2,inp3):
	record = Entrez.elink(dbfrom=inp1, db=inp2, id=inp3, retmax = 99999)
	return Entrez.read(record)

@retry(Exception)
def urlopen_with_retry(url): 
	return opener.open(url, timeout=15)

@retry2(Exception)
def urlopen_with_retry2(url): 
	return opener.open(url, timeout=15)

#map between uniprot or gids
def uniprot_mapping(fromtype, totype, identifier):
	base = 'http://www.uniprot.org'
	tool = 'mapping'
	params = {'from':fromtype,
				'to':totype,
				'format':'tab',
				'query':identifier,
	}
	data = urllib.urlencode(params)
	url = base+'/'+tool+'?'+data
	response = urlopen_with_retry(url)
	dat = StringIO.StringIO(response.read())
	page = np.loadtxt(dat,dtype=str)
	try:
		page = page[:,1][1:]
	except IndexError:
		return []
	return list(page)

#search db	
def esearch_id(db,inp):
	record = esearch_retry(db,inp)
	return record["IdList"]

#LINK between db1 to db2	
def elink_ids(inp,db1,db2):
	ids = []
	chunks = [inp[x:x+120] for x in xrange(0, len(inp), 120)]
	for chunk in chunks:
		record = elink_retry(db1,db2,",".join(map(str,chunk)))
		try:
			for i in record[0]['LinkSetDb'][0]['Link']:
				ids.append(i['Id'])
		except IndexError: pass
	return list(set(ids) - set(inp))

#send requests for CIDS from AIDS	
def aid_cid(inp,activity):
	ids = []
	chunks = [inp[x:x+40] for x in xrange(0, len(inp), 40)]
	for chunk in chunks:
		time.sleep(0.1)
		ids += cid_pug(chunk,activity)
		ids = list(set(ids))
	return list(set(ids) - set(inp))

#retrieve CIDS from AIDS	
def cid_pug(inp,activity):
	ids = []
	url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/" + ",".join(map(str,inp)) + "/cids/XML?cids_type=" + activity
	try:
		response = urlopen_with_retry2(url)
		root = ET.fromstring(response.read())
		for i in root.iterchildren():
			for x in i.iterchildren():
				ids.append(int(x.text))
	except: pass
	return ids

#retrieve SMILES from CIDS
def POSTactivities(inp):
    url = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
    #post CIDS to get	
    postdata = """
        <PCT-Data>
        <PCT-Data_input>
        <PCT-InputData>
        <PCT-InputData_download>
        <PCT-Download>
        <PCT-Download_uids>
        <PCT-QueryUids>
        <PCT-QueryUids_ids>
        <PCT-ID-List>
        <PCT-ID-List_db>pccompound</PCT-ID-List_db>
        <PCT-ID-List_uids>
        %(query)s
        </PCT-ID-List_uids>
        </PCT-ID-List>
        </PCT-QueryUids_ids>
        </PCT-QueryUids>
        </PCT-Download_uids>
        <PCT-Download_format value="smiles"/>
        <PCT-Download_compression value="gzip"/>
        </PCT-Download>
        </PCT-InputData_download>
        </PCT-InputData>
        </PCT-Data_input>
        </PCT-Data>"""
	#check SMILES from CIDS response	
    check_temp = '''
        <PCT-Data>
        <PCT-Data_input>
        <PCT-InputData>
        <PCT-InputData_request>
        <PCT-Request>
        <PCT-Request_reqid>%(reqid)s</PCT-Request_reqid>
        <PCT-Request_type value="status"/>
        </PCT-Request>
        </PCT-InputData_request>
        </PCT-InputData>
        </PCT-Data_input>
        </PCT-Data>'''

    def getfile(page):
        smileslist = []
        try:
            link = re.search( "<PCT-Download-URL_url>(.*)</PCT-Download-URL_url>", page).group(1)
        except AttributeError: print page
        response = urllib2.urlopen(link)
        compressedFile = StringIO.StringIO(response.read())
        response = None
        decompressedFile = gzip.GzipFile(fileobj=compressedFile)
        compressedFile = None
        ofile = decompressedFile.read().splitlines()
        return [x.split('\t') for x in ofile]

    smiles = []
    chunks = [inp[x:x+300000] for x in xrange(0, len(inp), 300000)]
    for i,chunk in enumerate(reversed(chunks)):
        percent = (float(i)/float(len(chunks)))*100 + 1
        sys.stdout.write(" Get Smiles %3d%%\r" % percent)
        sys.stdout.flush()
        smis = ['<PCT-ID-List_uids_E>' + str(ch) + '</PCT-ID-List_uids_E>' for ch in chunk]
        postd = postdata % {"query": "".join(map(str,smis))}
        req = urllib2.Request(url, postd, headers={'Content-Type': 'application/xml'})
        ret = urlopen_with_retry(req).read()
        while '<PCT-Status value="queued"/>' in ret or '<PCT-Status value="running"/>' in ret:
            time.sleep(3)
            m = re.search( "<PCT-Waiting_reqid>(\d+)</PCT-Waiting_reqid>", ret)
            if m: ret = urllib2.urlopen(url, check_temp % {"reqid": m.group(1)}).read()
        smiles += getfile(ret)
    return smiles

def main(targets):
	for targ in targets:
		print ' Processing target: ' + str(targ)
		#if egids get uniprots
		if options.egids:
			uniprots = uniprot_mapping("P_ENTREZGENEID","ACC",targ)
			print ' Mapped Uniprots: ' + ' '.join(map(str,uniprots))
			gids = []
			for u in uniprots:
					gids += uniprot_mapping("ACC","P_ENTREZGENEID",u)
					gids += esearch_id("gene",u)
		#else convert uniprot to egids
		else: gids = uniprot_mapping("ACC","P_ENTREZGENEID",targ)
		gids = list(set(gids))
		print ' Mapped GIDs: ' + ' '.join(map(str,gids))
		pids = elink_ids(gids,"gene","protein")
		print ' Mapped PIDs: ' + ' '.join(map(str,pids))
		aids = []
		aids += elink_ids(gids,"gene","pcassay")
		aids += elink_ids(pids,"protein","pcassay")
		aids = list(set(aids))
		print ' Number of AIDs: ' + str(len(aids))
		if len(aids) == 0: return []
		if options.actives:
			cids = aid_cid(aids,'active')
			of = open(targ+'.actives.smi','w')
			smiles = POSTactivities(cids)
			for s in smiles: of.write(' '.join(map(str,s))+'\n')
			of.close()
		if options.inactives:
			of = open(targ+'.inactives.smi','w')
			smiles = POSTactivities(cids)
			for s in smiles: of.write(' '.join(map(str,s))+'\n')
			of.close()
		return cids

if __name__ == '__main__':	
	if not options.actives or options.inactives:
		print 'Must supply --actives or --inactives! Exiting...'
		sys.exit()
	introMessage()
	main(open(options.inf,'r').read().splitlines())