# step 3
# read the pdb file and the ali file
# put it into the align2d function and get the pap file
# the pap file is the ali sequence and the pdb's sequence's comparation
from modeller import *
import os
import sys
import Bio
import urllib  
import shutil    
def comp(argv) :  
  env = environ()
#read the index for each type of the sequence to get the name of the pdb file and sequence file   
  if argv[1] == "PEMA" :
     file = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEMA/compare_test", "r")
  elif argv[1] == "PECA" :
     file = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA/compare_test", "r")
  elif argv[1] == "PEER" :
     file = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEER/compare_test", "r")
  for line in file :
# correct the format of the file so that we can read correctly, and use align2D to get the compare file  
     aln = alignment(env)
     column = line.split('\t')
     original = column[0]
     pdb = column[1].strip()
     pname = pdb[:-1]
     if pdb.len() >= 5 :
        ptype = pdb[4]
     else :
        ptype ='A'
     pfilename = pname+'.pdb'
     ofilename = original + '.ali'
     mdl = model(env, file=pname, model_segment=('FIRST:'+ptype,'LAST:'+ptype))
     aln.append_model(mdl, align_codes=pdb, atom_files=pfilename)
     aln.append(file=ofilename, align_codes=original)
     aln.align2d()

     aln.write(file=original + '-' + pdb +'.pap', alignment_format='PAP')
     oldpath = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/" + original + '-' + pdb +'.pap'
# clean the file create at the middle part but useless in the final part     
     if argv[1] == "PEMA" :
        if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEMA"):
               os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEMA")
        newpath = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEMA/" + original + '-' + pdb +'.pap'
        os.rename(oldpath,newpath)
     elif argv[1] == "PECA" :
        if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPECA"):
               os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPECA")
        newpath = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPECA/" + original + '-' + pdb +'.pap'
        os.rename(oldpath,newpath)
     elif argv[1] =="PEER" :
        if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEER"):
               os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEER")
        newpath = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEER/" + original + '-' + pdb +'.pap'
        os.rename(oldpath,newpath)
     os.remove(pfilename)
     os.remove(ofilename)
comp(sys.argv)