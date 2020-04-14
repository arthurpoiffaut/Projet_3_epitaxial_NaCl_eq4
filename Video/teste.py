import cv2
import numpy as np
import glob
import os
 #%%
img_array = []
path = r"C:\Users\Admin\Documents\Projet 3\Projet_3_epitaxial_NaCl_eq4\Video\videoFrames"
folder = os.listdir(path)

#%%
for filename in folder:
    img = cv2.imread(path + "\\" + filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)
 
#%%
out = cv2.VideoWriter('project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()