##### Read-me file for WAVES project
- code
    - bash
        - mothur.workflow.sh
            - This is a bash file that runs the mothur workflow we are using. The settings are as follows
            - The sequences are aligned to the 16S rRNA gene fragment from the Silva database.
            - Sequences were preclustered to a difference of 2
            - OTUs were classified using the Green Genes dataset
            - Subsampled to 4421, which is how many sequences B1014 had.
        - collect.files.sh
            - This is to move output files of interest from mothur into the dataEdited file for easier access.

- dataRaw - All the raw reads are currently stored on BensLabHD. They are also found on the EC2 Volume I have running.

- metadata - This contains several files from when Emily originally collected the samples. I haven't looked through them much.

- notes.md - This file contains my notes on working through the mothur pipeline the first time for this project.
