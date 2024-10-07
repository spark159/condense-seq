# Condense-seq

## Overview

Condense-seq is a high-throughput technique designed to directly measure the biophysical properties of native nucleosomes at a genome-wide scale with single-nucleosome resolution [1]. In this method, native mononucleosomes are purified from cells and condensed in vitro using various condensing agents, including polyamines, polycations, and heterochromatin proteins. By sequencing the supernatant phase of the condensate and comparing it to the input, the condensation propensity, or “condensability,” of individual genomic nucleosomes can be determined. Our findings demonstrate that native nucleosome condensability illuminates the biophysical principles driving genomic organization [1].

## Workflow

```mermaid
graph LR;
    classDef file fill:#FFF,stroke:#333,stroke-width:2px
    classDef script fill:#f9f,stroke:#333,stroke-width:2px

    Read1("Read1.fastq"):::file
    Read2("Read2.fastq"):::file
    Bowtie2("Bowtie2"):::script
    AlignedReads("AlignedReads.bam"):::file
    titr_file("TitrationFile.csv"):::file

    Read1 & Read2 --- Bowtie2
    Bowtie2 --> AlignedReads

    subgraph Quality control of data
        direction LR

        lendist_script("lendist.py"):::script
        lendist_file("ReadlengthFile.txt"):::file
        occ_script("occupancy.py"):::script
        occ_file("NCPoccupancyFile_occ.cn"):::file
        motif_script("motif.py"):::script
        motif_file("MotifFile_motif.txt"):::file

        lendist_script --> lendist_file
        occ_script --> occ_file
        motif_script --> motif_file
    
    end

    subgraph Genomic binning analysis
        direction LR

        bincount_script("bincount.py"):::script
        bincount_file("BincountFile_bin.cn"):::file

        bincount_script --> bincount_file

    end

    subgraph Single nucleosome analysis
        direction LR

        NCPpeak_script("NCPpeak.py"):::script    
        NCPpeak_file("NCPpeakFile_peak.cn"):::file

        NCPpeak_script --> NCPpeak_file

        NCPcov_script("NCPcov.py"):::script
        NCPcov_file("NCPcoverageFile_Ncov.cn"):::file
        
        NCPcov_script --> NCPcov_file

    end

    subgraph Sliding window analysis
        direction LR

        binsig_script("Binsig.py"):::script
        binsig_file("BinnedSignalFile_Bsig.cn"):::file

        binsig_script --> binsig_file

    end

    AlignedReads --- lendist_script
    AlignedReads --- occ_script
    NCPpeak_file --- motif_script

    coverage_script("coverage.py"):::script
    coverage_file("CoverageFile_cov.cn"):::file

    AlignedReads --- coverage_script
    coverage_script --> coverage_file

    AlignedReads --- bincount_script
    AlignedReads --- NCPpeak_script
    NCPpeak_file & coverage_file --- NCPcov_script
    coverage_file --- binsig_script

    NCPnum_script("NCPnum.py"):::script
    NCPnum_file("NCPnumberFile_num.cn"):::file

    bincount_file --- NCPnum_script    
    NCPcov_file --- NCPnum_script
    binsig_file --- NCPnum_script
    titr_file --- NCPnum_script

    NCPnum_script --> NCPnum_file

    NCPscore_script("NCPscore.py"):::script
    NCPscore_file("ScoreFile_score.cn"):::file

    NCPnum_file --- NCPscore_script
    NCPscore_script --> NCPscore_file

    Chalf_script("get_Chalf.py"):::script
    Chalf_file("ChalfFiles_Chalf.cn"):::file

    NCPnum_file --- Chalf_script
    Chalf_script --> Chalf_file

```

## Usage
<details>
<summary> bincount.py </summary>

Binning reference genome and get aligned read counts for each bin

  ```
  python bincount.py AlignedReads.bam -x ref_genome -w bin_size -o out_fname
  ```

</details>

<details>
<summary> coverage.py </summary>

Reading SAM/BAM files to get read coverage along reference genome.

  ```
  python coverage.py AlignedReads.bam -x ref_genome --chr chromosome -o out_fname --skip
  ```

</details>

</details>

<details>
<summary> NCPpeak.py </summary>

Peak calling for each nucleosome positions

  ```
  python NCPpeak.py AlignedReads.bam -x ref_genome --chr chromosome -o out_fname --skip
  ```

</details>

</details>

<details>
<summary> NCPcov.py </summary>

Compute coverage area under each nucleosome peaks

  ```
  python NCPcov.py NCPpeakFile_peak.cn CoverageFile_cov.cn --chr chromosome -o out_fname
  ```

</details>

</details>

<details>
<summary> Binsig.py </summary>

Compute coverage area for each sliding window along genome

  ```
  python Binsig.py CoverageFile_cov.cn -x ref_genome --Bsize bin_size --Bstep Bin_step --chr chromosome -o out_fname
  ```

</details>

<details>
<summary> NCPnum.py </summary>

Using titration file, estimate molecular number of nucleosomes for each bin or peak

  ```
  python NCPnum.py BincountFile_bin.cn | NCPcoverageFile_Ncov.cn | BinnedSigFile_Bsig.cn -t TitrationFile.csv --tnum TitrationNumber --chr chromosome -o out_fname
  ```

</details>

</details>

<details>
<summary> NCPscore.py </summary>

Get condensability score, which is a negative log of molecular number ratio over input, for each genomic bin or peaks

  ```
  python NCPscore.py NCPnumFile_num.cn --inpu NCPnumFile_num.cn -o out_fname
  ```

</details>

<details>
<summary> get_Chalf.py </summary>

Compute condensation point (C 1/2) by fitting logistic curve to molecular number changes over titrations

  ```
  python get_Chalf.py NCPnumFile_num.cn
  ```

</details>

## Publication

[1] https://www.biorxiv.org/content/10.1101/2023.12.08.570828v1
