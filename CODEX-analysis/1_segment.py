import os
import numpy as np
import skimage.io as io

from deepcell.utils.plot_utils import create_rgb_image
from deepcell.utils.plot_utils import make_outline_overlay
from deepcell.applications import Mesmer
from deepcell_toolbox.metrics import Metrics
from PIL import Image
from tqdm import tqdm
import matplotlib.pyplot as plt
app = Mesmer()

### Set the input parameters
img_path = ''     # Path to the image
fuse_channel = [1,14,18]    # Index of membrane channels
output_dir = 'output'   # Output directory
os.makedirs(output_dir, exist_ok=True)

img_data = io.imread(img_path)
img_data = np.transpose(img_data, (1,2,0))
print(img_data.shape)
img_input = np.zeros((img_data.shape[0], img_data.shape[1], 2))
img_input[:,:,0] = img_data[:,:,0]
max_fuse = np.max(img_data[:,:,fuse_channel], axis=2)
img_input[:,:,1] = max_fuse

X_train = img_input.reshape(1, img_input.shape[0], img_input.shape[1], img_input.shape[2])
print(X_train.shape)

segmentation_predictions = app.predict(X_train, image_mpp=0.5)
im = Image.fromarray(segmentation_predictions[0,:,:,0])
im.save(os.path.join(output_dir, 'cell_mesmer.tif'))

rgb_images = create_rgb_image(img_input.reshape(1, img_input.shape[0], img_input.shape[1], img_input.shape[2]), channel_colors=['blue', 'green'])
overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions)
fig, ax = plt.subplots(1, 2, figsize=(15, 15))
ax[0].imshow(rgb_images[0, ...])
ax[1].imshow(overlay_data[0, ...])
ax[0].set_title('Raw data')
ax[1].set_title('Cellular Predictions')
fig.savefig(os.path.join(output_dir,'_segmentation_cell_overlay.png'))