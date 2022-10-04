---
title: Setup
---
# Software Installation

Installing the software may take up to 30 minutes. You may also need to contact your local Information Technology Help Desk to get permission or assistance installing the software. Please do this **before the workshop**. We will not delay the start of the course while you install software. We *will* help you in advance to make sure that you have everything that you need.

Download and install the following software:

1. [R/4.2.1](https://cran.r-project.org/): Select the installation for your operating system (Windows, Mac, or Linux).
1. [RStudio](https://www.rstudio.com/products/rstudio/download/): Download the free **Rstudio Desktop**. 

Once you have installed R and RStudio, open RStudio to verify that the installation was successful.

# R Library Installation

In RStudio, copy and paste the following commands into the Console:

~~~
install.packages(c("tidyverse", "Seurat"), dependencies = TRUE)
~~~
{: .r}

Once the installation has finished, copy and paste the following commands into the console ot verify that both packages installed correctly.

~~~
library(tidyverse)
library(Seurat)
~~~
{: .r}

## Project Setup

1. Create a new project in your Desktop called `mapping`. 
    - Click the `File` menu button, then `New Project`.
    - Click `New Directory`. 
    - Click `New Project`.
    - Type `mapping` as the directory name. Browse to your Desktop to create the project there.
    - Click the `Create Project` button.

2. Use the `Files` tab to create  a `data` folder to hold the data, a `scripts` folder to house your scripts, and a `results` folder to hold results. Alternatively, you can copy and paste the following commands into the R console for step 2 only. You still need to create a project with step 1.

~~~
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
~~~
{: .r}

# Data Download

**Before the workshop**, please download the following files:

TBD: Users will download a subset of the counts and metadata from somewhere (Box?). So far, 25% of the full data from liveratlas.org seems to work well.

> DAS: Add some code here to download the file into `data`.

{% include links.md %}
