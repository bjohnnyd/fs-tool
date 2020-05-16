+++
title="Input Data"
description="Describes the main input data for the tool, especially netMHCpan predictions"
+++

The required input is a file of binding predictions for the HLA class I alleles and proteome of interest produced by NetMHCpan. 
The developers of netMHCpan were unable to provide a compatible download so this step needs to be performed separately by the user
and cannot be performed from command line by `fs-tool`. 

There are two ways to run NetMHCpan, either as a standalone command line tool or by accessing it over the server. 

Both versions are available at:

[https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.0](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.0)

The directions on obtaining the standalone version (for academic use only) are available under the `Downloads` tab.

To generate the input file using the web server, there are details under the `Instruction` tab, but briefly go to the 
`Submission` tab:

- Input your proteome of interest either via the input box or browsing to locate the fasta file on your local disk (one letter amino acid code).
- Select your peptide length. `fs-tool` can work with 8, 9, 10 or 11-mers. In the work presented in our paper we used 9-mers.
- Select your HLA class I alleles of interest, this should be all HLA class I alleles present in your cohort. 
NetMHCpan takes up to a maximum of 20 alleles per run so If you need more than 20 alleles then you will need to generate predictions for 20 alleles at a time, then concatenate the resulting binding predictions together.
- The thresholds set in the netMHCpan run will be used in the `fs-tool`.
- Choose whether to define binders based on eluted ligand predictions (the default) or binding affinity predictions (for this check the box). 
In our paper we choose to use binding affinity predictions rather than eluted ligand predictions since eluted peptides are thought to be heavily biased towards the most strongly binding peptides.
- Leave the “sort by prediction score” and “save predictions to XLS file” boxes  unchecked.
- Submit the run

Copy and paste or download the resulting raw output and save this locally and note the path as this will be your input file.

An example input file can be seen [HERE](https://raw.githubusercontent.com/bjohnnyd/fs-tool/master/tests/netmhcpan/netmhcpan_wBA.txt).




