import os
import numpy as np
import matplotlib.pyplot as plt 

# prepare the subplots
rows = 2
cols = 2
plt.subplots(rows,cols)

# Plot the Original Image
img = np.fromfile("Image",np.uint8).reshape(512,512)
plt.subplot(rows,cols,1)
plt.title("Original Image")
plt.imshow(img,cmap="gray", vmin=0, vmax=255)

# get available k-compressed images
ks = []
for file in os.listdir("build/"):
    if file.startswith("Image_"):
        ks.append(int(file.split('_')[1]))
ks = sorted(ks)

# Plot all k-compressed images
i = 2
for k in ks:
    img = np.fromfile(f'build/Image_{k}',np.uint8).reshape(512,512)
    plt.subplot(rows,cols,i)
    plt.title(f"Compressed Image k={k}")
    plt.imshow(img,cmap="gray", vmin=0, vmax=255)
    i += 1
    if i>rows*cols:
        break

plt.tight_layout()
plt.savefig("build/plot_exercise1.pdf")
