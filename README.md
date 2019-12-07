# Clustering
Project for exploration of clustering approaches

## Background
This repository contains Matei's Rmd lecture on clustering.  

I also recommend that you download the Young/Old FlowRepository dataset.  To do so:

1.  Install the FlowRepositoryR package
1.  Using that package, retrieve the __FR-FCM-ZZGS__ dataset
1.  Look at the code in this git repository as suggestion/guidance for how to
QC and pre-gate the data

After that, have fun!  Try FlowSOM.  Think of your own ideas for clustering
and population identification.  Look for T cell subsets that differ significantly
between young and old subjects.

## The Young/Old Dataset
This dataset was create by Rochester Human Immunology Center, David H. Smith Center 
for Vaccine Biology & Immunology, Rochester, NY (USA).  The purpose was to use SWIFT's 
competitive clustering assignment method to measure the differences between PBMC 
sub-populations in Old/Young subjects.  SWIFT is a very nice clustering method, 
implemented in MATLAB, described in [here](http://www2.ece.rochester.edu/projects/siplab/Software/SWIFT.html).  
Unfortunately, since it is implemented in MATLAB, which is a commercially available 
system, you may not be able to run it yourself.  However, you may be able to compare
your results to the published SWIFT results.

## QC and Gating
There are three files containing my R code in this repo.

1. yo_utils.R
    * This file contains utility functions for reading, gating, etc.
1. survey_data.R
    * Contains code to examine the data for warts, ultimately to be used to eliminate
    data files corrupted by un-correctable artifact.
1. gate_data.R
    * This file applies gating methods worked out in the surveying step, and writes
    gated FCS files to be used as input to downstream clustering workflows.

Please feel free to use these as suggestions for your own work.  Below is a summary of
my thought processes as I proceeded on this work.

### Evolution of the code and the thinking behind it
As I've previously mentioned, I like to declare functions in their own separate file:
in this case, __yo_utils.R__.  I select the "Source on Save" option in the editor,
and save frequently.  This has the effect of detecting any syntax errors as I go.

I then started writing the script __survey_data.R__.  The idea here was to look for any glaring
staining or other technical artifacts that might be too severe to correct - for example
by normalization.  At this point, yo_utils.R only had the first function, get_sample().
The remaining ones were added later.

I had made the decision to use SSC-A, CD3 and LIVEDEAD as
initial gating parameters.  So the first bit of this script (down to line 46 or so) simply made some figures
that I could examine manually.  These were the figures written into the pic_base
directory called raw_*.png.  These I could look at to get a first idea of how the
data were going to behave.

I then had the idea to use FlowFP to perform a more
quantitative survey of lack of conformity of individual samples relative to the
aggregate of all.  Here is where I went back and added ff_list on line 26 to accumulate
all data into a flowSet.  However, my 16 GB machine ran out of memory when I ran it,
so I realized that I needed to subsample the data (see line 34), randomly selecting
10,000 events from each flowFrame.  Then lines 48-54 compute and display a QC figure,
which clearly indicates a problem with instance #15.
![alt text](qc_raw.png)

At this point it was time to look at the other parameters.  I decided that the way
I would do this would be to first gate T cells, and then fingerprint all of the remaining
parameters.  To do this i needed to write the gating code.  So, back to __yo_utils.R__,
I thought a bit, then decided that the gating strategy would be as follows:

1. __gate_clean()__.  Look for "stationary" acquisition.  The assumption here is that
it shouldn't matter at what point during an acquisition an event was detected.  Should
be the same everywhere.  I selected SSC-A and one parameter off of each laser, with
the idea that if there was a fluidics hiccup I'd see it, and if there was a laser
hiccup on any laser I'd see that too.
1. __gate_singlet()__.  I used a trick I learned from Derrick Jones that just uses FSC-W
and SSC-W.  It's sort of a short-cut from the usual.
1. __gate_live()__.  I looked for the LIVEDEAD- and CD3+ blob.  This took a little fiddling
to get it to work for all files - in particular I pre-gated to allow blob.boundary()
to be able to ignore some of the junk down low.
1. __gate_scat()__.  This was a cleanup gate, since I'd already used CD3 in the previous
step.  This was particularly helpful in a few cases (e.g. instances 73, 101, 106, ...)
where the live gate snagged a few too many events with a bit of expression on the
live/dead marker.

The loop from lines 60-78 in __survey_data.R__ applied the gating strategy and stored a sub-sampled
version of the result.  The remaining lines repeated the flowFP QC calculation,
this time _including all fluorescence parameters not used for gating_.

Here you can see that, in addition to instance #15, there are several that should 
be checked out (which I haven't done as of this writing!.)  For example, the instances
colored orange to red are of some concern.

Finally, I wrote a new script, __gate_data.R__.  This script repeated the gating
calculations (without re-writing the gating figures), and output gated FCS files
to be used downstream.  Note that I would advise that you  __censor__ any instances that do not
meet your QC criteria.  A convenient way to do this is to create a spreadsheet
that contains one row for each FCS file.  It has a column, labeled for example
"valid".  Each file that you like gets a value of 1 for example, and files that
you don't like get a value of 0 (it doesn't matter what these values are, but
your downstream code will skip over the invalid files).

Hope this gives you some insight into how I go about starting a project!






