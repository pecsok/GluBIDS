#####################################################
*GluBIDS*: A Preprocessing Pipeline for GluCEST Data
#####################################################

*****
About
*****

.. image:: https://raw.githubusercontent.com/PennLINC/aslprep/main/docs/_static/aslprepworkflow.png

*ASLPrep* is a Arterial Spin Labeling  (ASL) data
preprocessing  and Cerebral Blood Flow (CBF) computation pipeline
that is designed to provide an easily accessible,
state-of-the-art interface that is robust to variations in scan acquisition
protocols and that requires minimal user input, while providing easily
interpretable and comprehensive error and output reporting.

It performs basic processing steps (coregistration, normalization, unwarping,
noise component extraction, segmentation, skullstripping etc.),
CBF computation, denoising CBF, CBF partial volume correction,
and providing outputs that can be easily submitted to a variety of group level analyses,
including task-based or resting-state CBF, graph theory measures, surface or volume-based statistics, etc.

The *GluBIDS* pipeline uses a combination of tools from well-known software
packages, including FSL_, ANTs_, FreeSurfer_ and AFNI_.
This pipeline was designed to provide the best software implementation for each state of preprocessing,
and will be updated as newer and better neuroimaging software become available.

This tool allows you to easily do the following:
