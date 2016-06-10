Hildebrand16 Guide
===================

###Introduction
The goal of this guide is to familiarize users with the Hildebrand16 electron micrograph (EM) dataset and the tools used to serve it. The dataset is a collection of serial-section electron microscopy image volumes acquired at different resolutions that collectively encompass the anterior quarter of one 5.5 days post-fertilization larval zebrafish.

We generated this dataset to study the anatomy of the larval zebrafish brain---which is included in its entirety---and peripheral nervous system. Toward this goal, we reconstructed and annotated some features of the nervous system, namely myelinated neuronal processes, which are also available for browsing or analysis. Many non-neuronal tissues are captured in the dataset.

----------

###Software
All electron micrographs and reconstructions are hosted using the [Collaborative Annotation Toolkit for Massive Amounts of Image Data (CATMAID)](http://catmaid.org/). CATMAID is designed to aid in the annotation and sharing of image datasets. This guide serves as an initial reference to help users navigate the CATMAID instance we use to host the dataset. For a more thorough understanding of the software and its functionality, please visit the [CATMAID documentation page](http://catmaid.readthedocs.org/).
Note that using the Google Chrome browser for interacting with CATMAID is highly recommended.

----------

###Getting started

####Accessing
Clicking on either the “[View data](http://hildebrand16.neurodata.io/catmaid/?pid=3&zp=537540&yp=351910&xp=303051&tool=navigator&sid0=2&s0=4)” or the "[View data with reconstructions](http://hildebrand16.neurodata.io/catmaid/?pid=3&zp=537540&yp=351910.65&xp=303051.44999999995&tool=tracingtool&sid0=2&s0=4)" link on the main page will immediately transport you into CATMAID.
The colored dots visible in the "View data with reconstructions" option are positions where a reconstructed object such as a myelinated neuron intersects with the current transverse section.

| View data        | View data and reconstructions | 
|:----------------:|:-----------------------------:|
| ![alt text][Vd]  | ![alt text][Vdar]             |
[Vd]: https://github.com/davidhildebrand/hildebrand16/raw/master/images/View_data_small.png "View data"
[Vdar]: https://github.com/davidhildebrand/hildebrand16/raw/master/images/View_data_and_reconstructions_small.png "View data and reconstructions"

####Navigating
Across the top of the page in CATMAID, you will see a toolbar. The toolbar provides useful information about the current view and contains various tools with which to interact with the data.
 
#####Orienting
Displayed in the toolbar are the current position, section number, and zoom level.
The current position is defined by the center of the field of view (x, y; in nm from the top left corner) and the section number (z-index; each ~60 nm):
![alt text][Tbl] 

The zoom level for the field of view can be set to visualize the full extents of the transverse section at low resolution or a restricted field of view at high resolution.
![alt text][Tbz]

The approximate resolution associated with each zoom level is:
| Zoom level | Resolution              | 
|:----------:|:-----------------------:|
| 5          | 1.8 μm/px               |
| 4          | 0.9 μm/px               |
| 3          | 451.2 nm/px             |
| 2          | 225.6 nm/px             |
| 1          | 112.8 nm/px             |
| 0          | 56.4 nm/px              |
| -1         | 28.2 nm/px              |
| -2         | 14.1 nm/px <sup>†</sup> |
| -3         | 7.05 nm/px <sup>†</sup> |
<sup>†</sup>Artificially upsampled. The actual resolution of raw data is 18.8 nm/px within most brain regions and 56.4 nm/px outside the brain.
[Tbl]: https://github.com/davidhildebrand/hildebrand16/raw/master/images/Toolbar_location.png "Toolbar location"
[Tbz]: https://github.com/davidhildebrand/hildebrand16/raw/master/images/Toolbar_zoom.png "Toolbar zoom"

#####Panning

To pan around the current transverse section, click the *pan tool button* on the toolbar:
![alt text][Tbp] 
Holding down your *left mouse button* while dragging will pan your view.
Alternatively, if you have a *middle mouse button*, hold it down when any tool is selected and drag to pan your view.
[Tbp]: https://github.com/davidhildebrand/hildebrand16/raw/master/images/Toolbar_pan.png "Toolbar pan"

#####Viewing reconstructions
To view existing neuron reconstructions, click the *tracing tool button* on the toolbar:
![alt text][Tbt] 
Colored dots will appear on top of the EM data:
![alt text][Vdar]
Each dot is a single node indicating where the current section intersects with a directed [polyline](https://en.wikipedia.org/wiki/Polyline) (or treeline) annotation used to represent the morphology of a particular neuron.
***ZOOM EXAMPLE***
[Tbt]: https://github.com/davidhildebrand/hildebrand16/raw/master/images/Toolbar_tracing.png "Toolbar tracing"

#####Viewing tags (meta-information)
To view the tags associated with individual nodes, click the *tag tool button* on the toolbar or press the *7 key*:
![alt text][Tbtg] 
***TAG EXAMPLE***

[Tbtg]: https://github.com/davidhildebrand/hildebrand16/raw/master/images/Toolbar_tags.png "Toolbar tags"

#####Finding reconstructed neurons

#####Viewing reconstructions in 3D

#####Different resolutions

#####Help
Pressing the *F1 key* will bring up a CATMAID help window pane that reveals commands available with any given tool. Note that some of these tools will not be available without additional access. For example, annotating additional features is not publicly available. However, access for creating additional annotations can be requested.

----------

> **Note:**

> - StackEdit is accessible offline after the application has been loaded for the first time.
> - Your local documents are not shared between different browsers or computers.
> - Clearing your browser's data may **delete all your local documents!** Make sure your documents are synchronized with **Google Drive** or **Dropbox** (check out the [<i class="icon-refresh"></i> Synchronization](#synchronization) section).



----------
Last updated on 2016-06-07 by David Hildebrand<!--se_discussion_list:{"h41SbNlsqb3mtPdQeOIdtotf":{"selectionStart":6194,"type":"conflict","selectionEnd":6204,"discussionIndex":"h41SbNlsqb3mtPdQeOIdtotf"}}-->