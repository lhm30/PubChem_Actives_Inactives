import urllib2
import string
import time
import urllib
from lxml import etree as ET
import csv
import time
from urlparse import urlparse
import httplib, sys
import time
import math
from functools import wraps
import csv
import gzip
import StringIO
import gc
import os
from Queue import Queue
from threading import Thread
import subprocess
from collections import OrderedDict
import datetime

class Worker(Thread):
    """Thread executing tasks from a given tasks queue"""
    def __init__(self, tasks):
        Thread.__init__(self)
        self.tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            func, args, kargs = self.tasks.get()
            try: func(*args, **kargs)
            except Exception, e: print e
            self.tasks.task_done()

class ThreadPool:
    """Pool of threads consuming tasks from a queue"""
    def __init__(self, num_threads):
        self.tasks = Queue(num_threads)
        for _ in range(num_threads): Worker(self.tasks)

    def add_task(self, func, *args, **kargs):
        """Add a task to the queue"""
        self.tasks.put((func, args, kargs))

    def wait_completion(self):
        """Wait for completion of all the tasks in the queue"""
        self.tasks.join()

concurrent = 6
opener = urllib2.build_opener()
opener.addheaders = [('User-agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_8_4) AppleWebKit/536.30.1 (KHTML, like Gecko) Version/6.0.5 Safari/536.30.1')]
opener.addheaders = [('Connection', 'keep-alive')]
opener.addheaders = [('Accept', 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8')]
opener.addheaders = [('Accept-Language', 'en-us')]

"""
Retry decorator used to resend url requests upon urllib failure, e.g: "[Errno 65] No route
to host", "[Errno 54] Connection reset by peer", "[Errno 8] nodename nor servname provided
or not known" which are caused by overwhelmed NCBI servers.

Each time the decorated function throws and exception, the decorator will wait a period of
time then retry calling the function until the maximum number of tries is used up.
@retry employs an exponential back-off doubling each wait time, e.g. '5sec,10,sec,20sec'

If the decorated function fails on the last try, the exception will occur unhandled.
"""
def retry(ExceptionToCheck, tries=4, delay=10, backoff=3, logger=None):
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
                        file = open('logger.txt', 'a')
                        file.write(msg)
                        file.close()
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)
        return f_retry
    return deco_retry

#attach retry decorator for GET links
@retry(Exception)
def urlopen_with_retry(url):
    return opener.open(url)

#attach retry decorator for POST data
@retry(Exception)
def post_urlopen_with_retry(url, postdata):
    req = urllib2.Request(url, postdata, headers={'Content-Type': 'application/xml'})
    return urllib2.urlopen(req)


"""
getGeneID function converts from 'Uniprot' accession to 'GENE ID' using the ESearch NCBI
eutils. The ESearch on the GENE database returns the GENE ID(s) for the Uniprot query.
"""

pool = ThreadPool(concurrent)

def getGeneID(term):
    gidlist = []
    #send uniprot term to esearch
    response = urlopen_with_retry('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=' + term + '&retmode=xml&RetMax=200')
    root = ET.fromstring(response.read())
    #search the xml for the gids
    for information in root.findall('IdList'):
        for gids in information.findall('Id'):
            gid=(gids.text)
            gidlist.append(gid)
    #remove duplicates
    gidlist = list(set(gidlist))
    return gidlist

aid_queue = Queue(concurrent*2)

"""
work for each thread from getAssayIDsFromGeneIDs
"""
def processAssayIDsFromGeneIDs(url):
    global aidlist
    response = urlopen_with_retry(url)
    context = ET.iterparse(StringIO.StringIO(response.read()), events=('end',), tag='Id')
    #search the xml for aids
    for event, elem in context:
        aidlist.add(elem.text)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    context = response = None

"""
getAssayIDsFromGeneIDs converts from 'GENE ID' to 'AID' using the ELink NCBI eutils.
The ELink retrieves PubChem BioAssays which have the GENE ID as an annotated target.
"""
def getAssayIDsFromGeneIDs(gids):
    global pool
    global aidlist
    global aid_queue
    #make threads for processing the gids
    try:
        #send blocks of 100 to be converted
        chunks = [gids[x:x+120] for x in xrange(0, len(gids), 120)]
        for chunk in chunks:
            formatted_gids = ''
            if len(gids) > 1:
                for countt, gid in enumerate(chunk):
                    if countt == 0:
                        formatted_gids = gid
                    else:
                        formatted_gids = formatted_gids + ',' + gid
            else:
                formatted_gids = gids[0]
            #put the chunks of gids in a queue
            url = ('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=pcassay&id=%s&retmode=xml') %(formatted_gids)
            pool.add_task(processAssayIDsFromGeneIDs, url)
        pool.wait_completion()
    except KeyboardInterrupt:
        sys.exit(1)
    for gid in gids:
        #remove the gid if it is the scraped aid list
        if gid in aidlist:
            aidlist.remove(gid)
    #remove duplicates
    return aidlist


"""
work for each thread from getProteinIDsFromGeneIDs
"""
def processProteinIDsFromGeneIDs(url):
    global pidlist
    response = urlopen_with_retry(url)
    context = ET.iterparse(StringIO.StringIO(response.read()), events=('end',), tag='Id')
    for event, elem in context:
        pidlist.add(elem.text)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    context = response = None

"""
getProteinIDsFromGeneIDs converts 'GENE ID' to 'PROTEIN ID' using the ELink NCBI eutils.
The ELink retrieves PubChem BioAssays which have the PROTEIN ID as an annotated target.
"""
def getProteinIDsFromGeneIDs(gids):
    global pool
    global pidlist
    global pid_queue
    try:
        #send blocks of 100 to be converted
        chunks = [gids[x:x+120] for x in xrange(0, len(gids), 120)]
        for chunk in chunks:
            formatted_gids = ''
            if len(gids) > 1:
                for countt, gid in enumerate(chunk):
                    if countt == 0:
                        formatted_gids = gid
                    else:
                        formatted_gids = formatted_gids + ',' + gid
            else:
                formatted_gids = gids[0]
            #send chunks of gid to elink to get proteins
            url = ('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id=%s&retmode=xml') %(formatted_gids)
            pool.add_task(processProteinIDsFromGeneIDs, url)
        pool.wait_completion()
    except KeyboardInterrupt:
        sys.exit(1)
    for gid in gids:
        #remove the gid if it is the scraped pid list
        if gid in pidlist:
            pidlist.remove(gid)
    #remove duplicates
    return pidlist

# <codecell>


"""
process the queue from getAssayIDsFromProteinIDs
"""
def processAssayIDsFromProteinIDs(url):
    global aidlist
    response = urlopen_with_retry(url)
    context = ET.iterparse(StringIO.StringIO(response.read()), events=('end',), tag='Id')
    for event, elem in context:
        aidlist.add(elem.text)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    context = response = None

"""
getAssayIDsFromProteinIDs converts from 'PROTEIN ID' to 'AID' using the ELink NCBI eutils.
The ELink retrieves PubChem BioAssays which have the PROTEIN ID as an annotated target.
"""
def getAssayIDsFromProteinIDs(pids):
    global pool
    global aidlist
    global aid_p_queue
    try:
        #send blocks of 100 to be converted
        chunks = [list(pids)[x:x+120] for x in xrange(0, len(pids), 120)]
        for chunk in chunks:
            formatted_pids = ''
            if len(pids) > 1:
                for countt, pid in enumerate(chunk):
                    if countt == 0:
                        formatted_pids = pid
                    else:
                        formatted_pids = formatted_pids + ',' + pid
            else:
                formatted_pids = pids[0]
            #send chunks of pids to elink to get assays
            url = ('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=pcassay&id=%s&retmode=xml') %(formatted_pids)
            pool.add_task(processAssayIDsFromProteinIDs, url)
        pool.wait_completion()
    except KeyboardInterrupt:
        sys.exit(1)
    for pid in pids:
        #remove the pid if it is the scraped aid list
        if pid in aidlist:
            aidlist.remove(pid)
    #remove duplicates
    aidlist = list(set(aidlist))
    return aidlist

# <codecell>

cidlist = set()

"""
process the queue from getInactiveCIDs
"""
def processInactiveCIDs(url):
    global cidlist
    response = urlopen_with_retry(url)
    context = ET.iterparse(StringIO.StringIO(response.read()), events=('end',), tag='Id')
    for event, elem in context:
        cidlist.add(elem.text)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    context = response = None

"""
getInactiveCIDs links from the comprehensive 'AID' list to the inactive PubChem
compounds (CID's) assayed in the studies.

The function does this using ELink and the linkname 'pcassay_pccompound_inactive' to
find the compounds 'tagged' as inactive for each AID study.
"""
def getInactiveCIDs(aids, uniprot):
    global uniprot_cid_list
    try:
        #send blocks of 40 to be converted
        #40 is optimised value. tests showed large xmls are slower due to size to parse
        chunks = [aids[x:x+60] for x in xrange(0, len(aids), 60)]
        for chunk in chunks:
            formatted_aids = ''
            if len(aids) > 1:
                for countt, aid in enumerate(chunk):
                    if countt == 0:
                        formatted_aids = aid
                    else:
                        formatted_aids = formatted_aids + ',' + aid
            else:
                formatted_aids = aids[0]
            #send chunks of aids to elink to get inactive compounds
            url = ('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pcassay&db=pccompound&id=%s&linkname=pcassay_pccompound_inactive') %(formatted_aids)
            pool.add_task(processInactiveCIDs, url)
        pool.wait_completion()
    except KeyboardInterrupt:
        sys.exit(1)
    for aid in aids:
        if aid in cidlist:
            cidlist.remove(aid)	
    for cid in cidlist:
        if not cid in uniprot_cid_list:
            uniprot_cid_list[cid] = []
        uniprot_cid_list[cid].append(uniprot)

"""
Once we have a list of the CIDs which are contained amoungst the final list of uniprots-
CID association. We post them using this function and then pass the ticket to the
getPOSTinactives function after this.
"""
def POSTinactives(cidterms):
    postdata = []
    ticket = ""
    postdata.append('<PCT-Data>')
    postdata.append('<PCT-Data_input>')
    postdata.append('<PCT-InputData>')
    postdata.append('<PCT-InputData_download>')
    postdata.append('<PCT-Download>')
    postdata.append('<PCT-Download_uids>')
    postdata.append('<PCT-QueryUids>')
    postdata.append('<PCT-QueryUids_ids>')
    postdata.append('<PCT-ID-List>')
    postdata.append('<PCT-ID-List_db>pccompound</PCT-ID-List_db>')
    postdata.append('<PCT-ID-List_uids>')
    for i in cidterms:
        postdata.append('<PCT-ID-List_uids_E>' + i + '</PCT-ID-List_uids_E>')
    postdata.append('</PCT-ID-List_uids>')
    postdata.append('</PCT-ID-List>')
    postdata.append('</PCT-QueryUids_ids>')
    postdata.append('</PCT-QueryUids>')
    postdata.append('</PCT-Download_uids>')
    postdata.append('<PCT-Download_format value="smiles"/>')
    postdata.append('<PCT-Download_compression value="gzip"/>')
    postdata.append('</PCT-Download>')
    postdata.append('</PCT-InputData_download>')
    postdata.append('</PCT-InputData>')
    postdata.append('</PCT-Data_input>')
    postdata.append('</PCT-Data>')
    postdata = ''.join(postdata)
    ticket = None
    response = post_urlopen_with_retry('http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi', postdata)
    response = response.read()
    context = ET.iterparse(StringIO.StringIO(response), events=('end',), tag='PCT-Waiting_reqid')
    for event, elem in context:
        ticket = elem.text
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
        break
    if ticket == None:
        context = ET.iterparse(StringIO.StringIO(response), events=('end',), tag='PCT-Download-URL_url')
        for event, elem in context:
            ticket = elem.text
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
            break
    if ticket == None:
        file = open('error_post.txt', 'a')
        file.write(postdata)
        file.close()
        file = open('error_return.txt', 'a')
        file.write(response)
        file.close()
        print "Couldn't find ticket. See error_post.txt & error_return.txt"
        sys.exit()
    context = response = None
    return ticket
# <codecell>

"""
We receive the ticket from the previous POSTinactives, and check is we have a response
from PUG. When we do we take the CID - smiles association.
"""
def getPOSTinactives(ticket):
    smileslist = []
    ok = False
    if not ticket.startswith("ftp://"):
        postdata = []
        link = ""
        postdata.append('<PCT-Data>')
        postdata.append('<PCT-Data_input>')
        postdata.append('<PCT-InputData>')
        postdata.append('<PCT-InputData_request>')
        postdata.append('<PCT-Request>')
        postdata.append('<PCT-Request_reqid>' + ticket + '</PCT-Request_reqid>')
        postdata.append('<PCT-Request_type value="status"/>')
        postdata.append('</PCT-Request>')
        postdata.append('</PCT-InputData_request>')
        postdata.append('</PCT-InputData>')
        postdata.append('</PCT-Data_input>')
        postdata.append('</PCT-Data>')
        postdata = ''.join(postdata)
        response = post_urlopen_with_retry('http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi', postdata)
        context = ET.iterparse(StringIO.StringIO(response.read()), events=('end',), tag='PCT-Download-URL_url')
        for event, elem in context:
            ok = True
            link = elem.text
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
            break
        context = response = None
        if ok == False:
            return None
    else :
        link = ticket
    response = urllib2.urlopen(link)
    compressedFile = StringIO.StringIO(response.read())
    response = None
    decompressedFile = gzip.GzipFile(fileobj=compressedFile)
    compressedFile = None
    file = decompressedFile.read().splitlines()
    for f in file:
        f = f.split("\t")
        smileslist.append(f[1])
    decompressedFile = None
    #print str(len(smileslist)) + " smiles dup"
    #we dont want to filter smileslist as will ruin the order and we want duplicate smiles
    #smileslist = list(set(smileslist))
    #smileslist = filter(bool, smileslist)
    #print str(len(smileslist)) + " smiles dup rem"
    return smileslist


aidlist = set()
pidlist = set()

uniprot_cid_list = OrderedDict()
original_uniprots = open('full_query_step.csv').read().splitlines()
entrycount = 1

#find inactive CIDs for uniprot query
print '******* Number of requested Uniprot IDs: ' + str(len(original_uniprots)) + ' *******\n'
for uniprot in original_uniprots:
    print 'Processing Entry No: ' + str(entrycount) + ' of ' + str(len(original_uniprots)) + ' - *** ' + uniprot + ' *** '
    gidlist = getGeneID(uniprot)
    aidlist = set()
    getAssayIDsFromGeneIDs(gidlist)
    pidlist = set()
    getProteinIDsFromGeneIDs(gidlist)
    gidlist = None
    getAssayIDsFromProteinIDs(pidlist)
    pidlist = None
    aidlist = list(set(aidlist))
    cidlist = set()
    getInactiveCIDs(aidlist, uniprot)
    entrycount += 1
    print str(len(cidlist)) + " inactives"
    file = open('inactivecount.txt', 'a')
    file.write(uniprot + "\t" + str(len(cidlist)) + "\n")
    file.close()

#send chunks of unique CIDs from the uniprot - CID association.

print "\n%s unique CIDs in dataset" %(str(len(uniprot_cid_list.keys())))
chunks = [uniprot_cid_list.keys()[x:x+500000] for x in xrange(0, len(uniprot_cid_list.keys()), 500000)]
tickets = {}
uniprot_smiles_list = {}
print "\nAquire SMILES"
for chunk_id, chunk in enumerate(chunks):
    tickets[chunk_id] = POSTinactives(chunk)
while len(tickets) > 0:
    for chunk_id, ticket in tickets.items():
        smileslist = getPOSTinactives(ticket)
        if smileslist != None:
            for index, smile in enumerate(smileslist):
                uniprot_smiles_list[smile] = uniprot_cid_list[chunks[chunk_id][index]]
                del uniprot_cid_list[chunks[chunk_id][index]]
            tickets.pop(chunk_id, None)
        if len(tickets) > 0:
			time.sleep(int(10/(len(tickets)+1)))

#Write smiles to file which will be passed ti cxcalc
file = open('smiles.txt', 'w')
for smile in uniprot_smiles_list.keys():
    file.write(smile + "\n")
file.close()

#determine the filter parameters for cxcalc
maxsize = 900
minsize = 125

"""
This function uses cxcalc mass to filter a list of smiles using preset thresholds.
"""
def filterSmiles(file):
	global uniprot_smiles_list
	original_size = len(uniprot_smiles_list)
	#remove large/small mass smiles	
	print '\nFiltering SMILES Size'
	#use the second cache to pass to 'cxcalc' and calculate the mass
	p2 = subprocess.Popen('cxcalc mass ' + file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#pipe the output to variable
	p2smiles = p2.stdout.readlines()
	file = open('filterdoutput.txt', 'a')
	for i in p2smiles:
		file.write(i + "\n")
	file.close()	
	#remove file headers		
	p2smiles.pop(0)
	listtodel = []
	#populate list of compounds too small/large to be deleted
	for s in p2smiles:
		#split line on tab to get [0] index [1] mass
		s = s.split('\t')
		#-1 from index to correct to position in python array
		rem = int(s[0])-1
		if float(s[1]) > maxsize:
			#if too large, tag the position in array to be deleted
			listtodel.append(rem)	
		if float(s[1]) < minsize:
			#if too small, tag the position in array to be deleted
			listtodel.append(rem)
		if float(s[1]) == 495.990:
			listtodel.append(rem)
			print "removed him!"
		if float(s[1]) == 542.082:
			listtodel.append(rem)
			print "removed him!"
	#remove the smiles to be deleted
	for delete in reversed(listtodel):
		#remove the smiles tagged for deletion from the previous array
		del uniprot_smiles_list[uniprot_smiles_list.keys()[delete]]
	final_size = len(uniprot_smiles_list)
	print 'Original size: ' + str(original_size) + '\tAfter deletion: ' + str(final_size) + '\tDeleted: ' + str(original_size - final_size) + '\tShould be: ' + str(len(listtodel)) 
	print '\nFiltered %s molecules\n' %(len(listtodel))
	return

"""
This function uses standardise mass to filter a list of smiles using preset thresholds.
"""	
def standardizeSmiles(file):

	print '\nStandardizing\n'	
	now = datetime.datetime.now()
	tfinish = len(uniprot_smiles_list)/12060
	then = now + datetime.timedelta(minutes = tfinish)
	print 'Estimated to finish at ' + str(then)
	#call the Chemaxon 'standardize', supply standardise settings, and pass the cache
	p1 = subprocess.Popen('standardize -c StandMoleProt.xml %s --log log.txt' %(file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#pipe the output to variable
	p1smiles = [x.strip() for x in p1.stdout.readlines()]
	
	return p1smiles

#filter the uniprot_smiles_list using the smiles we wrote to file
filterSmiles('smiles.txt')

#write filtered smiles to file
file = open('filtered_smiles.txt', 'w')
for smile in uniprot_smiles_list.keys():
    file.write(smile + "\n")
file.close()

stdsmiles = standardizeSmiles('filtered_smiles.txt')
print "\nStandardized " + str(len(stdsmiles)) + " SMILES"
file = open('stdizdoutput.txt', 'a')
for i in stdsmiles:
	file.write(i + "\n")
file.close()
	

#standardise the filtered smiles, write them to file and also update the cids in the
#uniprot smiles list to the final filtered, standardised smiles
file = open('standardized_smiles.txt', 'w')
for index, smile in enumerate(stdsmiles):
    for uniprot in uniprot_smiles_list.values()[index]:
		file.write(uniprot + "\t" + smile + "\n")
file.close()
		
print "\nFinished\n"