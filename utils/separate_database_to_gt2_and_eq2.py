import os
import tarfile
import shutil

'''
下の情報を編集後、databaseフォルダで実行すること！
python3 ./scripts/separate_database_to_gt2_and_eq2.py
'''

tmp_dir = os.path.abspath('./_tmp')
base_dir = os.path.abspath('.')
a_database_basename = 'pali2.8b'
a_database_ext = '.tar.gz'

# -------------------------------------
a_database_gt2 = os.path.join(base_dir,a_database_basename + '_gt2')
a_database_eq2 = os.path.join(base_dir,a_database_basename + '_eq2')

def unpack(a_archive):
	arc_file = tarfile.open(a_archive)
	arc_file.extractall(path=tmp_dir)
	arc_file.close()


try:
	os.makedirs(tmp_dir)
	os.makedirs(os.path.join(a_database_gt2,'inputdata'))
	os.makedirs(os.path.join(a_database_gt2,'ref'))
	os.makedirs(os.path.join(a_database_eq2,'inputdata'))
	os.makedirs(os.path.join(a_database_eq2,'ref'))
except:
	pass

unpack(os.path.join(base_dir, a_database_basename + a_database_ext))

for root,dirs,files in os.walk(os.path.join(tmp_dir, a_database_basename)):
	for a_file in files:
		if '.ref_fasta' in a_file:
			print('found {0}'.format(a_file))
			# count the number of sequences
			num_seq = 0
			with open(os.path.join(root,a_file)) as f:
				for a_line in f:
					if len(a_line) > 0 and a_line[0] == '>': num_seq = num_seq + 1

			# equal to 2
			if num_seq > 2:
				# *.ref_fasta
				shutil.move(os.path.join(root,a_file), os.path.join(a_database_gt2,'ref',a_file))
				# *.msf
				a_msf = a_file[:-len('.ref_fasta')] + '.msf'
				shutil.move(os.path.join(root,a_msf), os.path.join(a_database_gt2,'ref',a_msf))
				# ../inputdata/*.fasta
				a_fasta = a_file[:-len('.ref_fasta')] + '.fasta'
				shutil.move(os.path.join(root,'..','inputdata',a_fasta), os.path.join(a_database_gt2,'inputdata',a_fasta))
			else:
				# *.ref_fasta
				shutil.move(os.path.join(root,a_file), os.path.join(a_database_eq2,'ref',a_file))
				# *.msf
				a_msf = a_file[:-len('.ref_fasta')] + '.msf'
				shutil.move(os.path.join(root,a_msf), os.path.join(a_database_eq2,'ref',a_msf))
				# ../inputdata/*.fasta
				a_fasta = a_file[:-len('.ref_fasta')] + '.fasta'
				shutil.move(os.path.join(root,'..','inputdata',a_fasta), os.path.join(a_database_eq2,'inputdata',a_fasta))
