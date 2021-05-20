#!/usr/bin/python
import os
import shutil
import numpy as np
import threading
import subprocess
#my current directory
pwd = os.getcwd()
#launch makefile
os.chdir("../code/")
os.chdir(pwd)

#copy executable into subfolders
#shutil.copy("../code/main", "./gibbs/h0/.")
#shutil.copy("../code/main", "./gibbs/h02/.")
#shutil.copy("../code/main", "./metropolis/h0/.")
#shutil.copy("../code/main", "./metropolis/h02/.")


params_metro_h0={"nspin":50,
                 "j":1.,
                 "h":0.,
                 "metro":1,
                 "nblk":20,
                 "nstep":10000
}

params_metro_h02={"nspin":50,
                 "j":1.,
                 "h":0.02,
                 "metro":1,
                 "nblk":20,
                 "nstep":10000
}

params_gibbs_h0={"nspin":50,
                 "j":1.,
                 "h":0.,
                 "metro":0,
                 "nblk":20,
                 "nstep":10000
}
params_gibbs_h02={"nspin":50,
                 "j":1.,
                 "h":0.02,
                 "metro":0,
                 "nblk":20,
                 "nstep":10000
}

def job_launcher(pwd, path, temp_range, params,h, iterations):
    os.chdir(pwd+path)
    if h==0:
        f=open("../ene.dat","w")
        f.close()
        f=open("../heat.dat","w")
        f.close()
        f=open("../chi.dat","w")
        f.close()
    else:
        f=open("../mag.dat","w")
        f.close()
    for i in temp_range:
        f=open("input.dat", "w")
        f.write(str(i)+"\n")
        for key, value in params.items():
            f.write(str(value)+"\n")
        f.close()
        for j in range(iterations):
            if j == 0:
                subprocess.run(["./main", "--first"], stdout=subprocess.PIPE)
            else:
                subprocess.run(["./main", "--restart"], stdout=subprocess.PIPE)
        if h==0:
            #energy
            f=open("results/ene.dat", "r")
            lines=f.readlines()
            f.close()
            value=lines[-1].split()[-2]
            err=lines[-1].split()[-1]
            f=open("../ene.dat","a")
            f.write(str(i)+str(" ")+str(value)+str(" ")+str(err)+str("\n"))
            f.close()

            #heat
            f=open("results/heat.dat", "r")
            lines=f.readlines()
            f.close()
            value=lines[-1].split()[-2]
            err=lines[-1].split()[-1]
            f=open("../heat.dat","a")
            f.write(str(i)+str(" ")+str(value)+str(" ")+str(err)+str("\n"))
            f.close()


            #chi
            f=open("results/chi.dat", "r")
            lines=f.readlines()
            f.close()
            value=lines[-1].split()[-2]
            err=lines[-1].split()[-1]
            f=open("../chi.dat","a")
            f.write(str(i)+str(" ")+str(value)+str(" ")+str(err)+str("\n"))
            f.close()

        else:
            #mag
            f=open("results/mag.dat", "r")
            lines=f.readlines()
            f.close()
            value=lines[-1].split()[-2]
            err=lines[-1].split()[-1]
            f=open("../mag.dat","a")
            f.write(str(i)+str(" ")+str(value)+str(" ")+str(err)+str("\n"))
            f.close()



    os.chdir(pwd)

temp=np.linspace(0.5,2.,15)
pwd=os.getcwd()

print("Running with metropolis and h=0")
job_launcher(pwd,"/metropolis/0",temp, params_metro_h0, 0, 5)
print("Running with metropolis and h=0.02")
job_launcher(pwd, "/metropolis/0.02",temp, params_metro_h02, 1, 5)
print("Running with gibbs and h=0")
job_launcher(pwd, "/gibbs/0",temp, params_gibbs_h0, 0, 5)
print("Running with gibbs and h=0.02")
job_launcher(pwd, "/gibbs/0.02",temp, params_gibbs_h02, 1, 5)
