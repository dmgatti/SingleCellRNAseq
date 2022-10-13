---
title: Setup
---
# Software Installation

Installing the software may take up to 30 minutes. You may also need to contact 
your local Information Technology Help Desk to get permission or assistance 
installing the software. You may be able to install the applications using
JAX Self Service software.
Please do this **before the workshop**. We will not delay the start of the 
course while you install software. We *will* help you in advance to make sure 
that you have everything that you need.

If you do not already have `R` and `RStudio` installed, 
download and install the following software:

1. [R/4.2.1](https://cran.r-project.org/): Select the installation for your 
operating system (Windows, Mac, or Linux).
1. [RStudio](https://www.rstudio.com/products/rstudio/download/): Download the 
free **Rstudio Desktop**. 

You do not need to install this *exact* version of `R`, but it would be good to
make sure your `R` is relatively recent (say, updated within the past year).

Once you have installed `R` and `RStudio`, open `RStudio` to verify that the 
installation was successful.

# R Library Installation

In RStudio, copy and paste the following commands into the Console:

~~~
install.packages(c("tidyverse", "Seurat"), dependencies = TRUE)
~~~
{: .r}

Once the installation has finished, copy and paste the following commands into 
the console to verify that both packages installed correctly.

~~~
library(tidyverse)
library(Seurat)
~~~
{: .r}

## Project Setup

1. Create a new project called `scRNA`. 
    - Click the `File` menu button, then `New Project`.
    - Click `New Directory`. 
    - Click `New Project`.
    - Type `scRNA` as the directory name. Create the project anywhere you like,
      but don't forget where you put it!
    - Click the `Create Project` button.
    This will create a file called `scRNA.Rproj` in the directory you just 
    created. In the future you can double-click on this file to open 
    `RStudio` in this directory. This will be the easiest way to interact
    with the files/code you produce in this workshop.

2. Use the `Files` tab to create  a `data` folder to hold the data, a `scripts` 
folder to house your scripts, and a `results` folder to hold results. 
Alternatively, you can copy and paste the following commands into the `R` 
console for step 2 only. You still need to create a project with step 1.

~~~
dir.create("data")
dir.create("scripts")
dir.create("results")
~~~
{: .r}

# Data Download

**Before the workshop**, please download the following files:

TBD: Users will download a subset of the counts and metadata from somewhere (Box?). So far, 25% of the full data from liveratlas.org seems to work well. 

Open the `scRNA.Rproj` project.

~~~
download.file(url = 'https://thejacksonlaboratory.box.com/shared/static/vfe1bwyqtypxs6p5k4z0cw7z7jczyan1.zip', 
              destfile = 'data/mouseStSt_invivo.zip',
              method   = 'curl', 
              extra    = ' -L ')
download.file(url = 'https://thejacksonlaboratory.box.com/shared/static/b153ueu2lie3st760maj4zr9u0vp7o2t.zip', 
              destfile = 'data/mouseStSt_exvivo.zip',
              method   = 'curl', 
              extra    = ' -L ')
unzip(zipfile = 'data/mouseStSt_invivo.zip' ,
      exdir   = 'data')
unzip(zipfile = 'data/mouseStSt_exvivo.zip',
      exdir   = 'data')
~~~
{: .r}

{% include links.md %}
