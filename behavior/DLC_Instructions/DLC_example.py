# -*- coding: utf-8 -*-
"""# DLC Example
Welcome! This is an example DLC video analysis, by Lior de Marcas.

you can copy lines from here, and run them directly in ipython inside the anaconda prompt ("line-by-line") 

**Running "line-by-line" is recommended for unexperienced users.**

To simply run the entire file "as is", enter the conda DLC environment, navigate to where this file is, and type python + this file name (python DLC_example.py). 
You don't need to enter ipython! Make sure to prepare the network (see "Prepare network for analysing" in the manual) before running.

## Basics
### Entering DLC environment
You always need to enter the DLC environment before running anything. Simply open the Anaconda Prompt & run the next command (uncomment it first!)
"""

# conda activate DLC

"""If you are running the entire file (and not "line-by-line"), navigate to the folder containing this file, by running (uncommented):"""

# cd /d <path to folder containing this file>
# For example:
# cd /d D:\DLC_Instructions

"""If you are running "line-by-line", it is recommended to navigate to where you wish to save the results, by running (uncommented):"""

# cd /d <path to folder containing the results>
# For example:
# cd /d D:\DLC_Instructions

"""### Working with ipython
You need to run each command inside ipython. To enter ipython, simply uncomment & run:
"""

# ipython

"""To exit ipython:"""

# exit()

"""## Bureaucratic Stage
First, we need to write the path to the network & the video we want to analyse.
**Change the parameters here to match your computer & video!**
"""

net_path = r"D:\DLC_Instructions\Homecage_CFC_22labels_T1-LdM-2021-12-25\config.yaml"
# path to the config.yaml file of the network. We will use that as network identifier.

vid_path = r"D:\DLC_Instructions\lh49_CFC_Basler_Cam__22940788__20200327_101115219_fixedshort.mp4"
# path to the video you want to analyse. Any cropped and / or down-sampled version of this video will be placed in the same folder.

dest_path = r"D:\DLC_Instructions"
# path to the folder you wish to save the analysis results in.

downsample_height = 640
# height of the video after down-sampling. Width will rescale.

nMice = 1
# number of mice in your video. Our network is trained to deal with ip to 2 non-interacting mice.

"""In order to run it as a file, you need to change some controls here, to declare if you wish to check your video, down-sample and / or crop.
If you are running "line-by-line" (copying each line to ipython and running it), then ignore this step, and any "if" statement in this file.
"""

vid_check = True
# if you do not want to check for video corruption (see manual), change to "False"

Downsample = True
# if you do not want to down-sample the video at all, change to "False"

Crop = False
# if you want to crop, change to "True".
# Note that cropping require GUI that may bug a little - if it does you may get an error.

"""## Analysis Stage
The first thing we must do is import deeplabcut.
"""

import deeplabcut

"""### Check for corruption
Video corruption may occur as you move the video between computers.
We can check & fix videos that were corrupt.
"""

if vid_check:
  from deeplabcut.utils.auxfun_videos import VideoReader
  vid = VideoReader(vid_path)
  vid.check_integrity()

# you should not get any printed output from this.

#################################################
# if you do, you need to run the next steps (line-by-line, outside of ipython!)
# see manual for more info, and don't forget to uncomment!

# exit() # to exit ipython if you are still inside it
# cd /d <path_to_the_folder_holding_your_video>
# ffmpeg -i <video_name> -c:v h264 -crf 18 -preset fast <fixed_video_name>

# for example:
# exit() # to exit ipython if you are still inside it
# cd /d ??????
# ffmpeg -i ????? -c:v h264 -crf 18 -preset fast ?????

# Note that after you exit ipython, you will need to run the Bureaucratic step
# again, and change the vid_path to the fixed video's path.
#################################################

"""### Crop Video
If your mouses can only be in a small part of the video, crop any unneeded parts.
GUI used for cropping may bug, and prevent exit - 

In that case, abort using "ctrl + c" or by clicking the GUI "x", and then run this part again - it should work now.

If you are running as a file, this may cause an error - 
Simply run this part "line-by-line", and change the path of the video to the cropped one.
"""

if Crop:
  vid_path = deeplabcut.auxfun_videos.CropVideo(vid_path,useGUI = True)

"""### Down-sample Video
The smaller you video, the better the speed, with usually only small damage to accuracy.

It is recommended to down-sample to about 62.5% in each dimension - meaning that if your video was 1280 X 1024 (width X height), your new video will be 800 X 640.
"""

if Downsample:
  vid_path = deeplabcut.auxfun_videos.DownSampleVideo(vid_path,height = downsample_height)

"""### Analyse Video
First, make sure you changed the "config.yaml" & "inference_cfg.yaml" to match the number of mice to you video.
Our network wasn't trained on more then 2 non - interacting mice in the same video. See the manual, "Prepare network for analysing" for more info.

After you changed that, run:
"""

deeplabcut.analyze_videos(net_path,[vid_path], use_shelve=True, auto_track=True, destfolder = dest_path, identity_only=True, n_tracks = nMice)

"""To convert to csv, add an empty text file with the same name as the video to the folder with the results. Then run:"""

deeplabcut.analyze_videos_converth5_to_csv(dest_path, videotype = ".txt")