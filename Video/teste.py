import cv2
import numpy as np
import matplotlib.pyplot as plt
#%%
vidcap = cv2.VideoCapture('flocon.mp4')
print(vidcap.isOpened())
images = []
success = True
i=0
count = 0
while success == True:
    count +=1
    success,image = vidcap.read()
    if count % 10 ==0:
        images.append(image)
#%%
count = 0
for i in images:
    cv2.imwrite('videoFrames/allo%d.jpg' %count,i)
    count +=1
    
#%%
#for i in range(0,len(images)):
#    images[i] = cv2.cvtColor(images[i],cv2.COLOR_BGR2GRAY)[10:-10,250:1050]


image = images[64]
edges = cv2.Canny(image,800,800)

plt.figure(figsize = [10,14])
plt.imshow(edges,cmap = 'gray')
plt.title('Edge Image'), plt.xticks([]), plt.yticks([])
#plt.xlim([220,1050])
plt.show()