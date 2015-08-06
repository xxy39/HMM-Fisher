HMM-Fisher
==========
The HMM-Fisher program identifies differentially methylated (DM) CG sites (DMCs) and regions (DMRs) from whole genome and targeted bisulfite sequencing (BS) data. This approach first uses a hidden Markov chain to model the methylation signals to infer the methylation state as Not methylated (N), Partly methylated (P), and Fully methylated (F) at each CG site for each individual sample. Then the Fisher Exact test is used to identify differentially methylated CG sites. Third, identified DM CG sites are summarized into regions based on their status and distance. This program takes aligned BS data in multiple samples and outputs identified DM CG sites and regions.

HMM-Fisher requires R installed. Ideally it is run in a Linux/Unix system. This program includes one document and two folders:
_____________________________________________________________________________________________________________
HMM.Fisher.user.manual.pdf: A copy of the user manual

HMM.Fisher.manuscript.pdf: A copy of HMM-Fisehr manuscript

HMM.Fisher.code: A folder containing all R source code files used for HMM-Fisher.

Comparing_five_DM_methods.manuscript.pdf: A copy of the manuscript comparing five differential methylation detection methods

example.data: A folder containing all example input data, an example.script.txt for running HMM-Fisher, and the output files generated from the example.script.txt.
_____________________________________________________________________________________________________________
