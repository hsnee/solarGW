# solarGW
This branch is for preparing a tutorial in the form of a Jupyter/IPython notebook. 

This repository is for analysing LIGO data (we're currently analysing S6 data which is available publicly at losc.ligo.org) to investigate gravitational waves from our Sun.

Currently working on code in solarbayes/, a Bayesian approach in measuring the credibility of detecting gravitational waves from the Sun. Main routine at solarbayes/likelihood_job. 

Old code at simplistic/, a simplistic approach to detecting gravitational waves, from the cross-correlation of the signals in the two detecotors when pointing at the Sun. Main routine at simplistic/batchstrain_bg.py