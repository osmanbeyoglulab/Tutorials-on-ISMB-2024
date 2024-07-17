

# Virtual Tutorial 4: Computational Approaches for Identifying Context-Specific Transcription Factors using Single-Cell Multi-Omics Datasets (ISMB 2024)<!-- omit in toc -->

##   <a name='TableofContents--omitintoc--'></a>Table of Contents<!-- omit in toc -->

- [Overview](#overview)
    - [Basic principles behind TF activity inference methods](#basic-principles-behind-tf-activity-inference-methods)
    - [Overview of computational TF inference methods based on single cell omics](#overview-of-computational-tf-inference-methods-based-on-single-cell-omics)
- [Hands-on Tutorial](#hands-on-tutorial)
    - [Session 1 : Hands-on experience in applying tools and interpreting results using multiple TF activity inference methods using public scRNA-seq](#session-1--hands-on-experience-in-applying-tools-and-interpreting-results-using-multiple-tf-activity-inference-methods-using-public-scrna-seq)
    - [Session 2: Hands-on experience in applying tools and interpreting results using TF activity inference methods using public CITE-seq](#session-2-hands-on-experience-in-applying-tools-and-interpreting-results-using-tf-activity-inference-methods-using-public-cite-seq)
    - [Session 3: Hands-on experience in applying tools and interpreting results using multiple TF activity inference methods using public scATAC-seq and multiome](#session-3-hands-on-experience-in-applying-tools-and-interpreting-results-using-multiple-tf-activity-inference-methods-using-public-scatac-seq-and-multiome)
- [Environment Setup](#environment-setup)
  - [Installation of Conda](#installation-of-conda)
  - [Managing Environment](#managing-environment)
 
- [Introduction](#introduction)
- [Learning Objectives for Tutorial](#learning-objectives-for-tutorial)
- [Intended Audience and Level](#intended-audience-and-level)
- [Schedule](#schedule)
- [References](#References)

##  <a name='o=Overview'></a>Overview
####   <a name='Basic principles behind TF activity inference methods'></a><a href="https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024/blob/main/overview/ISMB 2024 Tutorial.pdf">Basic principles behind TF activity inference methods</a>

####   <a name='Overview of computational TF inference methods based on single cell omics'></a><a href="https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024/blob/main/overview/ISMB 2024 Tutorial.pdf">Overview of computational TF inference methods based on single cell omics</a>


##  <a name='Hands-onTutorial'></a>Hands-on Tutorial
####   <a name='Session1:Hands-onexperienceinapplyingtoolsandinterpretingresultsusingmultipleTFactivityinferencemethodsusingpublicscRNA-seq'></a><a href="https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024/blob/main/hands-on_tutorial/session-1/before_start.ipynb">Session 1 : Hands-on experience in applying tools and interpreting results using multiple TF activity inference methods using public scRNA-seq</a>


####  <a name='Session2:Hands-onexperienceinapplyingtoolsandinterpretingresultsusingTFactivityinferencemethodsusingpublicCITE-seq'></a><a href="https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024/blob/main/hands-on_tutorial/session-2/pySPaRTAN/pySPaRTAN_tutorial.ipynb">Session 2: Hands-on experience in applying tools and interpreting results using TF activity inference methods using public CITE-seq</a>

  

####  <a name='Session3:Hands-onexperienceinapplyingtoolsandinterpretingresultsusingmultipleTFactivityinferencemethodsusingpublicscATAC-seqandmultiome'></a><a href="https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024/blob/main/hands-on_tutorial/session-3/scenicplus.ipynb">Session 3: Hands-on experience in applying tools and interpreting results using multiple TF activity inference methods using public scATAC-seq and multiome</a>

##  <a name='EnvironmentSetup'></a>Environment Setup
###  <a name='InstallationofConda'></a>Installation of Conda

1. [Download the installer by choosing the proper installer for your machine.](https://www.anaconda.com/download/)
2. [Verify your installer hashes using SHA-256 checksums.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#hash-verification)
3. Install the installer:
	- Windows: Double-click the `.exe` file. Follow the instructions on the screen. For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html#installing-on-windows).
	- macOS: double-click the `.pkg` file. Follow the instructions on the screen. For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html#installing-on-macos).
	- Linux: In your terminal window, run: `bash Anaconda-latest-Linux-x86_64.sh`. Follow the prompts on the screen. For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#installing-on-linux).


###  <a name='ManagingEnvironment'></a>Managing Environment

With Conda, you can create, export, list, and update environments that have different versions of Python and/or packages installed in them. The JupyterLab, which can run in Conda environment,  is a web application for computational documents so that our code can produce rich, interactive output.


Below we will demonstrate how to create a Conda environment and install JupyterLab and packages for each tutorial session on macOS/Linux. You need to create a separate Conda environment for each session.

Use the **terminal** for the following steps. For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

#### Session-1:

1. Create an environment with python 3. 
   
        conda create --name stan python=3.12 -y
  
2. Activate the environment you just created: 
   
        conda activate stan
3. Install JupyterLab: 
   
        pip install jupyterlab

4. Install required packages

        git clone https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024.git

        cd Tutorials-on-ISMB-2024

        pip install -r requirements_stan.txt


#### Session-2:

1. Create an environment with python 3. 
   
        conda create --name pyspartan python=3.12 -y
  
2. Activate the environment you just created: 
   
        conda activate pyspartan
3. Install JupyterLab: 
   
        pip install jupyterlab

4. Install required packages

        git clone https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024.git

        cd Tutorials-on-ISMB-2024

        pip install -r requirements_pyspartan.txt


#### Session-3:

1. Create an environment with python 3. 
   
       conda create --name scenicplus python=3.11 -y
  
2. Activate the environment you just created: 
   
       conda activate scenicplus

3. Install JupyterLab: 
   
       pip install jupyterlab

4. Install required packages

       git clone https://github.com/aertslab/scenicplus

       cd scenicplus

       pip install .

       pip install dask==2024.5.0


After installing the required packages for the tutorial, launch JupyterLab and open the Jupyter notebook in each session subfolder to start the tutorial. To launch JupyterLab, enter the following command in your terminal:

      Jupyter lab


##  <a name='Background'></a>Introduction
Development of specialized cell types and their functions are controlled by external signals that initiate and propagate cell-type specific transcriptional programs. Activation or repression of genes by key combinations of transcription factors (TFs) drive these transcriptional programs and control cellular identity and functional state. For example, ectopic expression of the TF factors Oct4, Sox2, Klf4 and c-Myc are sufficient to reprogram fibroblasts into induced pluripotent stem cells. Conversely, disruption of TF activity can cause a broad range of diseases including cancer. Hence, identifying context-specific TFs is particularly relevant to human health and disease.

Systematically identifying key TFs for each cell-type represents a formidable challenge. Determination of TF activity in bulk tissue is confounded by cell-type heterogeneity. Single-cell technologies now measure different modalities from individual cells such as RNA, protein, and chromatin states. For example, recent technological breakthroughs have coupled the relatively sparse single cell RNA sequencing (scRNA-seq) signal with robust detection of highly abundant and well-characterized surface proteins using index sorting and barcoded antibodies such as cellular indexing of transcriptomes and epitopes by sequencing (CITE-seq). But these approaches are limited to surface proteins, whereas TFs are intracellular. Single-cell sequencing assay for transposase-accessible chromatin (scATAC-seq) measures genome-wide chromatin accessibility and reveals cellular memory and response to stimuli or developmental decisions. Recently several computational methods have leveraged these omics datasets to systematically estimate TF activity influencing cell states. We will cover these TF activity inference methods using scRNA-seq, scATAC-seq, Multiome and CITE-seq data through hybrid lectures and hand-on-training sessions. We will cover the principles underlying these methods, their assumptions and trade-offs. We will apply multiple methods, interpret results and discuss strategies for further in silico validation. The audience will be equipped with practical knowledge, essential skills to conduct TF activity inference independently on their own datasets and interpret results.
## <a name='LearningObjectivesforTutorial'></a>Learning Objectives for Tutorial
At the completion of the tutorial, participants will gain understanding into the basic concepts and recent advances in transcription factor inference methods for single-cell omics datasets including scRNA-seq, scATAC-seq, CITE-seq and Multiome. Four learning objectives are proposed:
1. Understand the basics principles underlying TF activity inference from single-cell omics
2. Understand the specific methodologies, assumptions, and trade-offs between computational inference methods
3. Gain hands-on experience in applying tools and interpreting results using multiple TF activity inference methods on public scRNA-seq, scATAC-seq, multiome and CITE-seq datasets
4. Discuss current bottlenecks, gaps in the field, and opportunities for future work.

## <a name='IntendedAudienceandLevel'></a>Intended Audience and Level
This tutorial is designed for individuals at the beginner to intermediate level, specifically targeting bioinformaticians or computational biologists with some prior experience in analyzing single-cell RNA sequencing (scRNA-seq), single-cell assay for transposase-accessible chromatin sequencing (scATAC-seq), Cellular Indexing of Transcriptomes and Epitopes by Sequencing (CITE-seq), and Multiome data, or those familiar with next-generation sequencing (NGS) methods. A foundational understanding of basic statistics is assumed.

While participants are expected to be beginners, a minimum level of experience in handling NGS datasets is required. The workshop will be conducted using Python and JupyterLab, necessitating prior proficiency in Python programming and familiarity with command-line tools.

To facilitate the learning process, participants will be provided with pre-processed count matrices derived from real datasets. All analyses, including JupyterLab notebooks and tutorial steps, will be available on GitHub for reference.

The tutorial will employ publicly accessible data, with examples showcased using datasets that will be made available through repositories such as the Gene Expression Omnibus or similar public platforms. This hands-on workshop aims to equip participants with practical skills and knowledge, enabling them to navigate and analyze complex datasets in the field of single-cell omics.


## <a name='Schedule'></a>Schedule
Tuesday, July 9, 2024 14:00 – 18:00 EDT
Time  | Tutorial
-------|-------------------
`14:00` | Welcome remarks and tutorial overview  <br /> 
`14:05` | Basic principles behind TF activity inference methods <br>  &nbsp; &nbsp; &nbsp;  * Overview of the importance of context-specific TF regulation in biological systems. <br> &nbsp;  &nbsp; &nbsp;  * Significance of TF dynamics in health and disease.<br> &nbsp;  &nbsp; &nbsp; * Single-cell multi-omics technologies for TF activity inference (scRNA-seq, scATAC-seq, Multiome<br /> 
`14:45` | Overview of computational TF inference methods based on single cell omics <br /> 
`15:45` | Break
`16:00` | Hands-on experience in applying tools and interpreting results using multiple TF activity inference methods using public scRNA-seq <br />
`16:45` | Hands-on experience in applying tools and interpreting results using TF activity inference methods using public CITE-seq <br />
`17:10` | Hands-on experience in applying tools and interpreting results using multiple TF activity inference methods using public scATAC-seq and multiome <br />
`17:55` |Discuss current bottlenecks, gaps in the field, and opportunities for future work <br />

##  <a name='References'></a>References
* Linan Zhang, April Sagan, Bin Qin, Baoli Hu, Hatice Ulku Osmanbeyoglu. STAN, a computational framework for inferring spatially informed transcription factor activity across cellular contexts: bioRxiv 2024.06.26.600782; doi: https://doi.org/10.1101/2024.06.26.600782

* Xiaojun Ma, Ashwin Somasundaram, Zengbiao Qi, Douglas J Hartman, Harinder Singh, Hatice Ulku Osmanbeyoglu, SPaRTAN, a computational framework for linking cell-surface receptors to transcriptional regulators, Nucleic Acids Research, Volume 49, Issue 17, 27 September 2021, Pages 9633–9647, https://doi.org/10.1093/nar/gkab745

* Kim, D., Tran, A., Kim, H.J. et al. Gene regulatory network reconstruction: harnessing the power of single-cell multi-omic data. npj Syst Biol Appl 9, 51 (2023). https://doi.org/10.1038/s41540-023-00312-6
