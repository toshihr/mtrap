# -*- coding: utf-8 -*-

'''

DATABASE FORMAT
./ref/*.ref_fasta  reference alignments
./in/*.fasta input sequence sets

USAGE EXAMPLE:
$> python generate_database.py prefab4 ../original_data/prefab4.tar.gz

'''

import os
import sys
import shlex
import tarfile
import subprocess
import shutil
import re
from xml.etree.ElementTree import *

try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO

def makedirs(dirName):
	try:
		os.makedirs(dirName)
	except OSError:
		pass

def split_text(text, n):
	#if not isinstance(text, unicode):
	#	text = text.decode('utf-8')
	io = StringIO(text)
	while True:
		s = io.read(n)
		if s:
			yield(s)
		else:
			break

def bali3_to_fasta(inFile, outFile):
	blocks = []

	tree = parse(inFile)
	elem = tree.getroot()

	coreblock = []
	for a_seq in elem.iter('colsco-data'):
		coreblock = a_seq.text.strip().split(' ')

	for a_name,a_original in zip(elem.iter('nid'),elem.iter('seq-data')):
		a_seq = ''
		for a_flag,a_char in zip(coreblock, a_original.text.strip()):
			if a_flag == '1':
				a_seq += a_char.upper()
			else:
				a_seq += a_char.lower()
		blocks.append( (a_name.text.strip(), a_seq) )

	with open(outFile, 'w') as fout:
		for a_name,a_seq in blocks:
			fout.write('>{0}\n'.format(a_name))
			for a_line in split_text(a_seq, 70):
				fout.write('{0}\n'.format(a_line))


def homstrad_to_fasta(inFile, outFile):
	blocks = []
	num = 0

	with open(inFile, 'rU') as fin:
		seq_buffer = []
		name = ''
		for line in fin:
			# skip line
			if line is None or len(line) == 0: continue
			if any([ len(line) >= len(a_str) and a_str == line[:len(a_str)] for a_str in ('structure','C;')]): continue

			# read
			if line[0] == '>':
				# store previous block
				if name != '':
					blocks.append( (name, ''.join(seq_buffer).replace('*','').replace('/','-') ) )
				# start new block
				name = line[1:].replace(';','_').strip()
				num += 1
				seq_buffer = []
			else:
				seq_buffer.append(line.strip().upper())

		if name != '':
			blocks.append( (name, ''.join(seq_buffer).replace('*','').replace('/','-') ) )

	with open(outFile, 'w') as fout:
		for a_name,a_seq in blocks:
			fout.write('>{0}\n'.format(a_name))
			for a_line in split_text(a_seq, 70):
				fout.write('{0}\n'.format(a_line))


def inFile_to_fastaFile(inFile, outFile):
	args = shlex.split('seqret -auto - supper1 -sequence {0} -outseq {1}'.format(inFile,outFile))
	try:
		p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		print('Failed to execute command: ' + args[0])
	# read outputs
	(stdoutdata,stderrdata) = p.communicate()

def fasta_to_ungappedFasta(inFile, outFile):
	blocks = []
	num = 0

	with open(inFile, 'rU') as fin:
		seq_buffer = []
		name = ''
		for line in fin:
			# skip line
			if line is None or len(line) == 0: continue

			# read
			if line[0] == '>':
				# store previous block
				if name != '':
					blocks.append( (name, ''.join(seq_buffer).replace('-','') ) )
				# start new block
				name = line[1:].strip()
				num += 1
				seq_buffer = []
			else:
				seq_buffer.append(line.strip().upper())

		if name != '':
			blocks.append( (name, ''.join(seq_buffer).replace('-','') ) )

	with open(outFile, 'w') as fout:
		for a_name,a_seq in blocks:
			fout.write('>{0}\n'.format(a_name))
			for a_line in split_text(a_seq, 70):
				fout.write('{0}\n'.format(a_line))


def unpack(arcName, outDir):
	print('unzip ' + arcName + ' to ' + outDir)
	# make a tmp directory
	makedirs(outDir)

	arc_file = tarfile.open(arcName)

	for t in arc_file:
		if t.isdir():
			# 755: rwxr-xr-x
			t.mode = 0o755
		else:
			# 644: rw-r--r--
			t.mode = 0o644

	arc_file.extractall(path=outDir)
	arc_file.close()

def pack(dirName, erase=False, addDirName=True):
	'''
		pack dir.
		ex.) ./[dirName]/data1
                             data2
                             ...
               are stored in dirName.tar.bz2 with relative path dirName/data1,...

	'''
	# .tar.gz -> mode=gz
	arcName = dirName + '.tar.bz2'
	th = tarfile.open(arcName, 'w:bz2')
	for root, dirs, files in os.walk(dirName):
		for f in files:
			fullpath = os.path.join(root, f)
			#above_dir = os.path.basename(root)
			#relativepath = os.path.join(above_dir, f)
			if addDirName:
				relativepath = fullpath[len(os.path.dirname(dirName) + os.path.sep):]
			else:
				relativepath = fullpath[len(dirName + os.path.sep):]
				
			# archive relative path
			th.add(fullpath, arcname=relativepath)
	th.close()

	if erase:
		shutil.rmtree(dirName, ignore_errors=True)

	return arcName

def make_database(arcName, arcType, arcSubname, convertFunc, outDir, tmpDir, refDir, refRegex, refExt, addRelativePath):
	arcName = os.path.abspath(arcName)

	# --- check ---
	if arcSubname not in arcName:
		print('arc name is wrong.')
		sys.exit(1)

	# --- unpack original database ---
	originalDir = os.path.join(tmpDir,arcType)
	unpack(arcName=arcName, outDir=originalDir)

	# --- original database structure ---
	refDir = originalDir if refDir is None else os.path.join(originalDir,refDir)

	# --- new database structure ---
	outRefDir = os.path.join(outDir,'ref')
	outInDir = os.path.join(outDir,'in')

	# --- reformat ---
	for root, dirs, files in os.walk(refDir):
		for f in files:
			# --- ext check ---
			if not ((refExt is None) or (refExt and len(f) > len(refExt) and f[-len(refExt):] == refExt)): continue

			# --- regex check ---
			if refRegex is not None:
				m = re.search(refRegex, f)
				if m is None: continue

			print('target: {0}'.format(f))

			fullpath = os.path.join(root, f)
			relativepath = fullpath[len(refDir + os.path.sep):]

			basename,ext = os.path.splitext(f)

			if not addRelativePath: relativepath = ''

			outRefFasta = os.path.join(outRefDir, os.path.dirname(relativepath), basename + '.ref_fasta')
			outInFasta = os.path.join(outInDir, os.path.dirname(relativepath), basename + '.fasta')

			# --- convert ---
			exec('{0}(inFile=fullpath, outFile=outRefFasta)'.format(convertFunc))

			fasta_to_ungappedFasta(inFile=outRefFasta, outFile=outInFasta)

	shutil.rmtree(originalDir)


if __name__ == '__main__':

	arcType_information = {
		'prefab1': dict(arcSubname='prefab1', refDir='refalns',refRegex=None , refExt=None, convertFunc='inFile_to_fastaFile', addRelativePath=False),
		'prefab2': dict(arcSubname='prefab2', refDir='refaln',refRegex=None ,refExt=None, convertFunc='inFile_to_fastaFile',addRelativePath=False),
		'prefab3': dict(arcSubname='prefab3', refDir='ref',refRegex=None ,refExt=None, convertFunc='inFile_to_fastaFile',addRelativePath=False),
		'prefab4': dict(arcSubname='prefab4', refDir='ref',refRegex=None ,refExt=None, convertFunc='inFile_to_fastaFile',addRelativePath=False),
		'homstrad': dict(arcSubname='homstrad_ali_only', refDir=None,refRegex=None ,refExt='.ali', convertFunc='homstrad_to_fasta', addRelativePath=False),
		'balibase3rv11': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV11'),refRegex='^BB[^S].*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv11s': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV11'),refRegex='^BBS.*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv12': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV12'),refRegex='^BB[^S].*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv12s': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV12'),refRegex='^BBS.*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv20': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV20'),refRegex='^BB[^S].*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv20s': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV20'),refRegex='^BBS.*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv30': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV30'),refRegex='^BB[^S].*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv30s': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV30'),refRegex='^BBS.*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv40': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV40'),refRegex='^BB[^S].*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv50': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV50'),refRegex='^BB[^S].*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'balibase3rv50s': dict(arcSubname='BAliBASE3', refDir=os.path.join('bb3_release','RV50'),refRegex='^BBS.*' ,refExt='.xml', convertFunc='bali3_to_fasta', addRelativePath=False),
		'sabmark_sup': dict(arcSubname='SABmark', refDir=os.path.join('SABmark','sup'),refRegex='^[^g][^r][^o][^u][^p]' ,refExt='.fasta', convertFunc='inFile_to_fastaFile', addRelativePath=False),
		'sabmark_twi': dict(arcSubname='SABmark', refDir=os.path.join('SABmark','twi'),refRegex='^[^g][^r][^o][^u][^p]' ,refExt='.fasta', convertFunc='inFile_to_fastaFile', addRelativePath=False),
	}

	if len(sys.argv) <= 2:
		print('generate_database [target type] [target arc name]')
		print('support target type: prefab[1-4], homstrad, balibase3rv[11,11s,12,12s,20,20s,30,30s,40,50,50s')
		sys.exit(1)

	targetType = sys.argv[1]
	targetArcName = sys.argv[2]

	outArcName = os.path.join('.',targetType + '.tar.bz2')
	if os.path.exists(outArcName):
		os.remove(outArcName)
		#print('there already exists {0} !'.format(outArcName))
		#sys.exit(1)

	tmpDir = os.path.join(os.environ['HOME'],'tmp')
	outDir = os.path.join(tmpDir, 'output')

	# clean
	if os.path.exists(outDir):
		shutil.rmtree(outDir)

	# make dir
	makedirs(os.path.join(outDir,'ref'))
	makedirs(os.path.join(outDir,'in'))

	if targetType in list(arcType_information.keys()):
		make_database(arcName=targetArcName, arcType=targetType, outDir=outDir, tmpDir=tmpDir, **arcType_information[targetType])
		os.rename(pack(dirName=outDir,erase=False,addDirName=False), outArcName)
		print('generate {0}'.format(outArcName))

	shutil.rmtree(outDir)