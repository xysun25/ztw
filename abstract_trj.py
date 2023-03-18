from operator import ne
from re import A
import sys
import os
from tokenize import Number
import numpy as np
import math
import matplotlib.pyplot as plt
from sympy import numbered_symbols

class read_trj():
    def __init__(self, trj_file):
        self.trj_file=trj_file


    def read_trj(self):
        frames=0
        timesteps=[]
  
        if os.path.exists(self.trj_file):
            f2=open(self.trj_file,'r')
            sum_of_atoms=f2.readline()
            sum_of_atoms=int(sum_of_atoms)
            for line in f2:
                if "Atoms." in line:
                    frames+=1
                    timesteps.append(line)
            f2.close()
            return frames,timesteps
        else:
            raise IOError("the trjactory file is not existed!")

    def abstract_trj(self,frames):
        sum_of_atoms=0
        timesteps=[]
        f3=open(self.trj_file,'r')
        sum_of_atoms=f3.readline()
        sum_of_atoms=int(sum_of_atoms)
        f3.seek(0)
        f4=open("abstracted_traj",'a+')
        frames=[int(x) for x in frames]
        print(frames)
        f,timesteps=self.read_trj()

        for line in f3:
            for i in frames:
                timestep=timesteps[i-1]
                #print(timestep)
                if timestep in line:
                    f4.write(str(sum_of_atoms)+'\n')
                    f4.write(timestep)
                    for i in range(sum_of_atoms):
                        line=f3.readline()
                        f4.write(line)

#read molecular information from mmol files
if os.path.exists('abstracted_traj'):
    os.remove("./abstracted_traj")
trj_file=input('please input the trj_file name:')
read_trj_inst=read_trj(trj_file)
frames,steps=read_trj_inst.read_trj()
frames=int(frames)
print(frames)

abstracted_frames=input("please input the frames you want to choose: all, one or range such as all or 1 2 3, or 1-5: ")
if abstracted_frames=="all":
    new_list=np.arange(1,frames+1,1)
    new_list.tolist()
    print(new_list)
    read_trj_inst.abstract_trj(new_list)
    
elif "-" in abstracted_frames:
    list=abstracted_frames.split("-")
    new_list=np.arange(int(list[0]),int(list[1])+1,1)
    new_list=new_list.tolist()
    print(new_list)
    read_trj_inst.abstract_trj(new_list)

else:
    new_list=abstracted_frames.split(' ')
    print(new_list)
    read_trj_inst.abstract_trj(new_list)



