# step 1
# this script is used to split the original ali file
# it divides the input file into 3 folders so that it can be easy for the next step
# 
import Bio.SeqIO
import os
import argparse
import sys
# to make the input file spirit into one sequence one file
def divide_file(argv):
# use the argument to select which kind of animals that we can get the sequence, for example PEMA
    pname = argv[1]
    output1 = ""
    output2 = ""
    path = ""
    ali = []
    name = []
    file = []
    i = 0
    if pname == "PEMA" :


      with open("/home/unhTW/share/mcbs913_2020/omics/scripts/PEMA_modller_input.ali", "r") as fullpema:
          for line in fullpema:
              i += 1

              if i % 3 == 1:
                  column = line.split(";")
                  nowname = column[1].strip()
                  nowname += '.ali'
                  name += column[1]
                  file += '>P1;'+column[1]

              else:
                  file += [line]
              if i % 3 == 0:
                  output1 = "".join(file)
          #create the path to store the split file for each sequence        
                  path = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali1",nowname)

                  if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali1"):
                     os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali1")
                  split = open(path.format(path), "w")
                  split.write(output1)
                  split.close()
                  path =""
                  output1 = ""
                  
                  file.clear()
                   
   
      fullpema.close()
      print("part 1\n")
      output2 = "".join(name)    

      try:
          with open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali1/name.txt".format(path), "w") as subset:
                  subset.write(output2)

      except IOError:
          print("IO error2 makeing new file")
          # this part is for the PECA
    elif pname =="PECA" :
      with open("/home/unhTW/share/mcbs913_2020/omics/scripts/PECA_modller_input.ali", "r") as fullpema:
          for line in fullpema:
              i += 1

              if i % 3 == 1:
                  column = line.split(";")
                  nowname = column[1].strip()
                  nowname += '.ali'
                  name += column[1]
                  file += '>P1;'+column[1]

              else:
                  file += [line]
              if i % 3 == 0:
                  output1 = "".join(file)
                  path = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali2",nowname)

                  if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali2"):
                     os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali2")
                  split = open(path.format(path), "w")
                  split.write(output1)
                  split.close()
                  path =""
                  output1 = ""
                  nowname = ""
                  file.clear()
                   
   
      fullpema.close()
      print("part 2\n")
      output2 = "".join(name)    
#create a index so that we can read it in a easy way
      try:
          with open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali2/name.txt".format(path), "w") as subset:
                  subset.write(output2)

      except IOError:
          print("IO error2 makeing new file")
          # this part is for the PEER
    elif pname == "PEER" :
      with open("/home/unhTW/share/mcbs913_2020/omics/scripts/PEER_modller_input.ali", "r") as fullpema:
          for line in fullpema:
              i += 1
#             print(i)
              if i % 3 == 1:
                  column = line.split(";")
                  nowname = column[1].strip()
                  nowname += '.ali'
                  name += column[1]
                  file += '>P1;'+column[1]
#             print(i)
              else:
                  file += [line]
              if i % 3 == 0:
                  output1 = "".join(file)
                  path = os.path.join("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali3",nowname)
#                 print(output1)
                  if not os.path.exists("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali3"):
                     os.makedirs("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali3")
                  split = open(path.format(path), "w")
                  split.write(output1)
                  split.close()
                  path =""
                  output1 = ""
                  
                  file.clear()
             #  if i == 3:
#                   break     
   
      fullpema.close()
      print("part 3\n")
      output2 = "".join(name)    

      try:
          with open("/home/unhTW/share/mcbs913_2020/omics/test_files_for_modeller/split_ali3/name.txt".format(path), "w") as subset:
                  subset.write(output2)

      except IOError:
          print("IO error2 makeing new file")
    else:
      name = "hello"
      name += ".ali"
      print("please input PEMA, PECA, PEER as argv[1] ")
#       split = open(name, "w")
#     except IOError:
#        print("IO error3 with gff file")


divide_file(sys.argv)
