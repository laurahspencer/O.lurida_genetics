### Preparing microsatellite data for analysis in GenePop 

#### Checking 2016/2017 Fidalgo Bay raw data for correct binning

Crystal rounded the raw microsat data for the Fidalgo Bay 2016-hatchery and 2017-wild data.  She provided both raw and rounded data. Before moving forward with the rounded data, I'll check out the binning method she used. 

In the Excel file [Olympic Oyster NFH_NFW (1).xlsx](https://github.com/laurahspencer/O.lurida_genetics/blob/master/Data/Olympic%20Oyster%20NFH_NFW%20(1).xlsx) she includes data from both wild and hatchery NF samples. She houses raw data for each locus in separate tabs, creates a list of "bins" at 0.2 increments, calculates frequencies for each bin and visualizes with histograms. 

![image](https://user-images.githubusercontent.com/17264765/35253966-0be0dc42-ff9d-11e7-85bb-15158f6f4335.png)

Then, using the frequency distributions she assigned alleles, for example: 

![image](https://user-images.githubusercontent.com/17264765/35254029-57188304-ff9d-11e7-812e-6a6268acd746.png)   |    ![image](https://user-images.githubusercontent.com/17264765/35254078-99f39bbe-ff9d-11e7-8541-74505b89351b.png)

One question I have is regarding the assignment of all even-numbered alleles for Oly10, Oly11 & Oly12, while alleles are odd for Oly13, Oly15 & Oly19. 

**I also noticed that Oly18 data was initially processed, then not completed nor included in the "rounded" tab.** I emailed Crystal to see what's up (I presume it was an oversight). 

Next step is to export the data into a GenePop format.  GenePop is one of the most commonly used programs used to analyze microsatellite data. There are several ways to use GenePop: [on the web, at the command line](http://genepop.curtin.edu.au/), and in [R](https://cran.r-project.org/web/packages/genepop/index.html).  I like to work in R. I could not find an R-based function to convert .csv format to GenePop format, however thre is an Excel plug-in caled [GenAlEx](http://biology-assets.anu.edu.au/GenAlEx/Download.html) that one can use.  I download version 6.503 (Dec 5, 2016). Then, I merged the wild and hatchery data into one spreadsheet. I also found online that the commonly used "genind" format has a few key formatting requirements, which I point out in the following screenshot: 

![snip20180122_24](https://user-images.githubusercontent.com/17264765/35255706-3575b87c-ffa5-11e7-8729-e23c08b4e870.png)

With this merged file open, I also opened the GenAlEx program. Then, I used the GenAlEx plug-in to export the file as a GenPop formatted .txt file: 

![image](https://user-images.githubusercontent.com/17264765/35255731-512f0dca-ffa5-11e7-8212-c378a2252891.png)

A window pops up, which should automatically ID the #loci & #samples if you formatted the spreadsheet like I did; I edited the "Title" to include 2016/2017 info.

![image](https://user-images.githubusercontent.com/17264765/35255755-7000ef70-ffa5-11e7-9edc-6062c98284f0.png)

Saved the file as a .txt file under [Oly2016NFH+2017NFW_Merged.txt](https://github.com/laurahspencer/O.lurida_genetics/blob/master/Data/Oly2016NFH%2B2017NFW_Merged.txt); here's what the resulting PopGen formatted file looks like: 

![image](https://user-images.githubusercontent.com/17264765/35255846-eb99b96e-ffa5-11e7-9364-e96df4e2cbc0.png)
