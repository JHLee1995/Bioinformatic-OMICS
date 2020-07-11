from modeller import *
import os
import sys
import Bio
import time
import urllib
import shutil
from Bio.PDB import PDBList
from io import StringIO

start_time = time.time()

def buildprofile(argv):
  log.verbose()
  env = environ()
  filename = StringIO()
  i = 0
  filecontent =[]
#   filename.write(argv[1])
#-- Prepare the input files

#-- Read in the sequence database
  sdb = sequence_db(env)
  sdb.read(seq_database_file='/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/basic-example/pdb_95.pir', seq_database_format='PIR',
           chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

#-- Write the sequence database in binary form
  sdb.write(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
            chains_list='ALL')

#-- Now, read in the binary database
  sdb.read(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
           chains_list='ALL')
  if argv[1] == "PEMA" :
    subset = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali1/name.txt", "r")
    print("step 1")
  elif argv[1] == "PECA" :
    subset = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali2/name.txt", "r")
  elif argv[1] == "PEER" :
    subset = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali3/name.txt", "r")
  for line in subset:
      pname = line.strip()
      name = pname + '.ali'
      i += 1
#-- Read in the target sequence/alignment
      if argv[1] == "PEMA" :
         path = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali1",name)
      elif argv[1] == "PECA" :
         path = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali2",name)
      elif argv[1] == "PEER" :
         path = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali3",name)
      aln = alignment(env)
      aln.append(file=path, alignment_format='PIR', align_codes='ALL',allow_alternates = True)

      prf = aln.to_profile()
      prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
                gap_penalties_1d=(-500, -50), n_prof_iterations=1,
                check_profile=False, max_aln_evalue=0)

      aln = prf.to_alignment()
      j = 0
      maxname = ""
      for seq in aln:
         tvl = aln[0]
         identity = seq.get_sequence_identity(tvl)
         # print(seq.atom_file,identity)
         if j == 1:
            max = identity
         elif j > 1 :
            if max < identity :
               max = identity
               maxname = seq.atom_file
         j += 1
      print(maxname)
      if maxname != "" :

         '''Selecting structures from PDB'''
         if len(maxname) > 4 :
            filename = maxname[:-1]
            type = maxname[4]
            print(filename,type)
         pdbl = PDBList()
         pdbl.retrieve_pdb_file(filename,pdir='PDB_files',file_format ='pdb')
         bpdbname = 'pdb'+filename +'.ent'
         pdbname = filename +'.pdb'
         bpdbpath = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/PDB_files",bpdbname)
         pdbpath = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller",pdbname)
         newpath = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller",name)
         shutil.copy(path,newpath)
         shutil.copy(bpdbpath,pdbpath)
         filecontent +=[pname + '\t' + maxname+'\n']

         # mdl = model(env, file=filename, model_segment=('FIRST:A','LAST:A'))
#          aln.append_model(mdl, align_codes=maxname, atom_files=pdbpath)
#          aln.append(file=path, align_codes=pname)
#          print(path)
#          aln.align2d()
#          aln.write(file=pname+'-'+filename+'.ali', alignment_format='PIR')
#          aln.write(file=pname+'-'+filename+'.pap', alignment_format='PAP')
 #   urllib.urlcleanup()

      print(i)
      if i == 20:
         break
  print("here")
  if argv[1] == "PEMA" :
     if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEMA"):
            os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEMA")
     output = "".join(filecontent)
     print(output)
     compare = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEMA/compare_test", "w")
     compare.write(output)
  elif argv[1] == "PECA" :
     if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA"):
            os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA")
     output = "".join(filecontent)
     compare = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA/compare_test", "w")
     compare.write(output)
  elif argv[1] == "PEER" :
     if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEER"):
            os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEER")
     output = "".join(filecontent)
     compare = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA/compare_test", "w")
     compare.write(output)
buildprofile(sys.argv)
print("--- %s seconds ---" % (time.time() - start_time))
