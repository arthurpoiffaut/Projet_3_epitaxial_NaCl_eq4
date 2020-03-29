import cv2
import numpy as np
import matplotlib.pyplot as plt
#%%
NumberOfFrames = 678-1 # frames
MaxRadius = 435.5633307157709 # pixels
FinalTime = 44*60 # secondes
FinalLength = 2.5*10**-3 #metres

frame2seconds = FinalTime/NumberOfFrames
pixel2meters = FinalLength/MaxRadius

#%%
vidcap = cv2.VideoCapture('flocon.mp4')
print(vidcap.isOpened())
images = []
frameVec = []
c = []
success = True
i=0
count = 0
step = 2
while success == True:
    count +=1
    success,image = vidcap.read()
    if count % step ==0:
        frameVec.append(count)
        images.append(image)
#%%        
frameVec = np.array(frameVec[0:-1])  

gray = [None]*(len(images)-1)
edges = [None]*(len(images)-1)
for i in range(0,len(images)-1):
    gray[i] = cv2.cvtColor(images[i],cv2.COLOR_BGR2GRAY)[10:-10,250:1050]
    edges[i] = cv2.Canny(gray[i],100,100)
    
 
#%%
count = 0
for i in edges:
    cv2.imwrite('videoFrames/allo%d.jpg' %count,i)
    count +=1
#%%
def findCenter(image):
    xMat = np.zeros(np.shape(image))
    yMat = np.zeros(np.shape(image))
    for i in range(0,np.shape(image)[0]):
        for j in range(0,np.shape(image)[1]):
            if image[i,j] == 0:
                xMat[i,j] = np.nan
                yMat[i,j] = np.nan
            else:            
                xMat[i,j] = i*image[i,j]
                yMat[i,j] = j*image[i,j]           
    center = np.array([np.nanmean(yMat),np.nanmean(xMat)]) 
    return center

def findRadius(image,center):   
    distance = np.zeros(np.shape(image))
    for i in range(0,np.shape(image)[0]):
        for j in range(0,np.shape(image)[1]):
            if image[i,j] == 0:
                distance[i,j] = np.nan
            else:
                distance[i,j] = np.sqrt((i-center[0])**2 + (j-center[1])**2)           
    radius = np.nanmax(distance)     
    return radius    


radius = [None]*len(edges)
center = [None]*len(edges)
growth = [None]*len(edges)
for i in range(0,len(edges)):
    print("step {0:} of {1:}".format(i,len(edges)))
    center[i] = findCenter(edges[i]/255)
    radius[i] = findRadius(edges[i]/255,center[i])
    if i == 0:
        growth[i] = 0
    else:
        growth[i] = radius[i] - radius[i-1]
growth = np.array(growth) 
#%%
#index = 10
#
#image = edges[index]
#
#plt.figure(figsize = [10,14])
#plt.imshow(image,cmap = 'gray')
#plt.scatter(center[index][0],center[index][1],color = 'red')
#plt.title('Edge Image')#, plt.xticks([]), plt.yticks([])
##plt.xlim([220,1050])
#plt.show()

plt.figure(figsize = [14,10])
plt.scatter(frameVec*(frame2seconds) , growth*(pixel2meters/frame2seconds)/step*10**6)    
plt.show()

avgGrowthRate = np.nanmean(growth*(pixel2meters/frame2seconds)/step)
print('The average growth rate is {:.2f} µm/s'.format(avgGrowthRate*10**6))

#%% write edges video
filename = "edges.avi"
fourcc = cv2.VideoWriter_fourcc(*'xvid')#no idea what this is
fps = 15
size = np.shape(edges[0])
video = cv2.VideoWriter(filename, fourcc, fps, size)
for i in range(0,len(edges)):
    video.write(images[i])
video.release()