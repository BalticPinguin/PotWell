#!/usr/bin/python3.4
import numpy as np
import sys, mmap, re

def main(argv):
    #read the file and store in mmap-object
    infile=open(argv[0], 'r')
    f=mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)
    infile.close()
    coordx=re.findall(b'(?<=coordx =)[,\-e\d .\n+]+', f)[0]
    xvec=re.findall(b"[\d.\-\+e]+",coordx)
    coordy=re.findall(b'(?<=coordy =)[,\-e\d .\n+]+', f)[0]
    yvec=re.findall(b"[\d.\-\+e]+",coordy)
    #yvec=coordy.split()
    #index=np.sort(yvec)
    z=False
    if re.search(b'(?<=coordz =)[,\-e\d .\n]', f) is not None:
        coordz=re.findall(b'(?<=coordz =)[,\-e\d .\n+]+', f)[0]
        zvec=re.findall(b"[\d.\-\+e]+",coordz)
        z=True
        # sort: y,z,x-changes...
        index=np.lexsort((zvec,xvec,yvec)) # ascending sorting f
    else:
        index=np.lexsort((xvec,yvec)) # ascending sorting f
    # now, I have the coordinates of the points.
    #  elem_num_map / node_num_map  > one of them connects the elements to coordinates??
    #HERE I WILL ASSUME UNITY-MAPPING FOR elem_num_map (as is in the shown file)
    #  vals_nod_var1 =  maps values to the elements!?
    numVar=re.findall(b'(?<=vals_nod_var)[\d ]+', f)
    Variables=[]
    #print(len(numVar))
    #print("")
    #print("")
    for i in range(len(numVar)//2): # because it appears in the header as well
        vari=re.findall(b'(?<=vals_nod_var)[=\,\-e\d .\n+]+', f)[i+len(numVar)//2] 
        #print(vari)
        Variables.append(re.findall(b"[e\d.\-\+]+",vari)[1:])
        #print("")
        #print("")
        #print(Variables)
        #print("")
        #print("")
    #print(np.shape(Variables))
    for i in range(len(xvec)):
        if i>0 and yvec[index[i]]!=yvec[index[i-1]]:
            print('')
        print("%8.5g"%(float(xvec[index[i]])),end="  ")
        print("%8.5g"%(float(yvec[index[i]])),end="  ")
        if z:
            print("%8.5g"%(float(zvec[index[i]])),end="  ")
        for j in range(len(numVar)//2):
            print("%8.5g"%(float(Variables[j][index[i]])),end="  ")
        print('')

if __name__ == '__main__':
   main(sys.argv[1:])
