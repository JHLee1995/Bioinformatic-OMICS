# step 4
# use the Orthogroups.csv to find the vcf file linked to the ali sequence
# check the sequence in the pap and know the target sequence is existed or not
# if exist, link it to the pdb file and get the atom place then find out the protein number
# then get the 3d factor
from modeller import *
import os
import sys
import Bio
import urllib  
import shutil  
# use the argument to select the type of the sequence that we need to use
def vcf_to_pdb(argv) :
   RA = []
   pdb = []
   match =""
   if argv[1] == "PEMA" :
      pdbfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEMA/compare_test", "r")
   elif argv[1] == "PECA" :
      pdbfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA/compare_test", "r")
   elif argv[1] == "PEER" :
      pdbfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEER/compare_test", "r")
# find the target pdb and the target vcf file name   
   for line in pdbfile :
      column = line.split('\t')
      RA += column[0]
      ranow = column[0]
      raid = ranow[-9:]
      
      pdb += column[1]
      pdbnow = column[1]
# read this file as an index to get the target vcf file and set the file name     
      index = open("/home/unhTW/share/mcbs913_2020/omics/files/protein_files/Results_Apr15/Orthogroups.csv","r")
      for content in index :
         col = content.split('\t',1)
         fname = col[0]
        
         Rnameset = col[1]
         
         Rname = Rnameset.strip()
         rname = Rname.split(',')
         
         for rid in rname :
            vcrid = rid[-9:]
            
            if vcrid == raid :
               str = fname + '\t' + ranow + '-' + pdbnow 
               
               match += str
      index.close()
   print(match)
# read the vcf and pdb compared file to get the file name   
   if argv[1] == "PEMA" :
      vtobfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEMA/vcf_to_pdb", "w")
      vtobfile.write(match)
      vtobfile.close()
      vtobrfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEMA/vcf_to_pdb", "r")
   elif argv[1] == "PECA" :
      vtobfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA/vcf_to_pdb", "w")
      vtobfile.write(match)
      vtobfile.close()
      vtobrfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPECA/vcf_to_pdb", "r")
   elif argv[1] == "PEER" :
      vtobfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEER/vcf_to_pdb", "w")
      vtobfile.write(match)
      vtobfile.close()
      vtobrfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/cPEER/vcf_to_pdb", "r")
# read the two file and check the protein is existed or not      
   for e in vtobrfile:
      s = e.split('\t')
      fname = s[0]+'.fa.vcf'

      papname = s[1].strip() + '.pap'
# read the pap file
      fpath = "/home/unhTW/share/mcbs913_2020/omics/scripts/vcf/"+fname
      if argv[1] == "PEMA" :
         ppath = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEMA/" + papname
      if argv[1] == "PECA" :
         ppath = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPECA/" + papname
      if argv[1] == "PEER" :
         ppath = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/papPEER/" + papname
      ffile = open(fpath,"r")
      pfile = open(ppath,"r")
      i = 0
      pos =[]
      for fline in ffile :
         i += 1
         if i > 4 :
            c = fline.split('\t')

            nownum = c[1]
            pos += [nownum]
      j = 0
      count = 0
      pdbsequence = []
      for pline in pfile :
         j += 1
         if j % 5 == 2 :
            p = pline.strip()
            pc = p[-45:]
            print(pc,"+=")
            pdbsequence += pc
      pdbid = s[1]
      pdbid = pdbid[-6:-2]
    
      oldpdb = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/PDB_files/pdb" + pdbid +".ent"
      newpdb = "/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/PDB_files/" + pdbid +".pdb"
      shutil.copy(oldpdb,newpdb)
      output =""
      pdbfile = open(newpdb,"r")
# get the target protein and find the target atoms 3d factor, so that the output file can be used into spanning tree code      
      for num in pos :
         no = int(num)
# if the protein is existed, read the pdb file.  
         if no < len(pdbsequence) :
            target = pdbsequence[no]
            
         else :
            target = '-'
         
         if target != '-' and target !=' ' :
#             print(target,"+++")
            for pdbline in pdbfile :
               pcol = pdbline.strip()
               
               if pcol[:4] =="ATOM" :
                  pcolumn = pcol.split(' ')
                  while '' in pcolumn :
                     pcolumn.remove('')
#                   print(pcolumn)
                  if pcolumn[5] == num :
                     print(pcolumn)
                     sline = pcolumn[5] + '\t' + pcolumn[6] +'\t' + pcolumn[7] + '\t' + pcolumn[8] +'\n'
                     
                     output += sline
      factorfile = open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/"+pdbid,"w")
      factorfile.write(output)
vcf_to_pdb(sys.argv)              