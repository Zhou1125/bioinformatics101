# bioinformatics101




Just released today 08th June 2021. An ITol tool alternative in R

- [Paper](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msab166/6294410) | [Usage](https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html) | [Github](https://github.com/YuLab-SMU/ggtreeExtra/)












2021 June 02 | [Gentle introduction on PacBio Hifi](https://youtu.be/1aqnLQ-nwYk)

MCScan code


I will use diamond for blastp and synvisio as visualization tool.

- Very important point 1: eliminate * symbole from the protein file


- Vey important point 2: Arrange the GGF file following the example from the github page in data


```
at1	AT1G01010	3631	5899

```

For A.majus use **mRNA** 

```
am8	Am08g37680.T01	3631	5899


```



fOR SESAMUM INDICUM use **CDS**

```
si1 SIN_1021768

```
For Lycopersicum use **mRNA**


```
sl3 Sol 
```




- Very important point 3: Rename the chromosome id by doinf AT1 for Arabidopsis thaliana chromosome 1


- Very important point 4: Performe the blastp run for every single pairing like if you have twop genomes A and B, do the blastp for A vs A, B vs B, A vs B and B vs A

- Very important point 5: for the GFF file eliminate with excel the duplicated genes 
- Very important point 6: Concatenate all blastp file as weell as gff files 
- After running use the output file for synvisio online. It is much more faster


Note: In some protein file, there is "*" . Diamond tools is sensitive to those star symboles. So It is important to remove them by doing:


```bash

cat final_perilla4x.protein_rename.fasta | sed 's/[*]//g ' >  final_perilla4x.protein_rename_cleaned.fasta


```



After getting the result of blast, you may need to convert the csv format into a tabulate format before running mcscan. So do:


```bash

sed 's|,|\t|g' pc_pf.blast.csv > pc_pf.blast

```




# Step 1: blastp  of pc versus pc

```python

source activate diamond_env

diamond makedb --in perilla_v1.0_protein_without_point.fasta  -p 64 -d perilla_v1.0_protein_without_point

diamond blastp -d perilla_v1.0_protein_without_point -q perilla_v1.0_protein_without_point.fasta -p 64 --evalue 0.00001 --out pc_vs_pc_diamond.cleaned.csv --outfmt 6 &> log.run.diamond.cleaned.pc_pc &

```

# Step 2: blastp of pc versus at


```python

diamond makedb --in Athaliana_447_Araport11.protein.cleaned.fa  -p 64 -d Athaliana_447_Araport11.protein.cleaned


diamond blastp -d Athaliana_447_Araport11.protein.cleaned -q perilla_v1.0_protein_without_point.fasta -p 64 --evalue 0.00001 --out pc_vs_ara_diamond.cleaned.csv --outfmt 6 &> log.run.diamond.cleaned.pc_ara &

```


# Step 3: blastp of at versus pc

```python

diamond makedb --in perilla_v1.0_protein_without_point.fasta  -p 64 -d perilla_v1.0_protein_without_point

diamond blastp -d perilla_v1.0_protein_without_point -q Athaliana_447_Araport11.protein.cleaned.fa -p 64 --evalue 0.00001 --out ara_vs_pc_diamond.cleaned.csv --outfmt 6 &> log.run.diamond.cleaned.ara_pc &

```

# step 4:blastp at versus at

```python
diamond makedb --in Athaliana_447_Araport11.protein.cleaned.fa  -p 64 -d Athaliana_447_Araport11.protein.cleaned


diamond blastp -d Athaliana_447_Araport11.protein.cleaned -q Athaliana_447_Araport11.protein.cleaned.fa -p 64 --evalue 0.00001 --out ara_vs_ara_diamond.cleaned.csv --outfmt 6 &> log.run.diamond.cleaned.ara_ara &

```



# Step 5: Concatenate both

```python

cat pc_vs_pc_diamond.cleaned.csv pc_vs_ara_diamond.cleaned.csv ara_vs_pc_diamond.cleaned.csv ara_vs_ara_diamond.cleaned.csv > pc_at.blast

```























Just see on twitter new tool for read mapping using long-reads | [Vulcan](https://gitlab.com/treangenlab/vulcan) | Use [vulcan](https://www.biorxiv.org/content/10.1101/2021.05.29.446291v1?s=03) for mapping and [sniffles](https://github.com/fritzsedlazeck/Sniffles/wiki) for SV calling |

Un bon [tutorial](https://ressources.france-bioinformatique.fr/sites/default/files/5%20-%20Variants%20Structuraux.pdf) sour les SV.





```

conda install -c bioconda vulcan

conda install -c bioconda sniffles


./vulcan -r ./test/GCF_000146045.2_R64_genomic.fna -i ./test/test_reads.fa -w ./test/ -o vulcan


./sniffles -m mapped.sort.bam -v output.vcf

```


Awesome disovring today 30 May 2021 | [Video](https://youtu.be/P18iAOj9LEk) | [bioarchive](https://www.biorxiv.org/content/10.1101/2021.03.12.435103v1.full) | [code avaliablity](https://github.com/UMNKimballLab/NWRGenomeAssembly_v1.0)


Each figure est bien expliquee avec les codes a l'appuie. Jai sauve les codes ici dans le github. Awesome.






TASSEL workshop  [video1](https://youtu.be/BY_6qSIW4Rk) | [video2](https://youtu.be/hiReCEINEVs) | [google group](https://groups.google.com/g/tassel) | [github](https://github.com/SushanRu/TASSEL-workshop) | [
University of Minnesota Plant Breeding Center](https://www.youtube.com/channel/UCl1JgjWhfTqG7DqRl6j_zlA)










From this [Assembling vertebrate genomes](https://www.nature.com/articles/s41576-021-00379-z) Arang Rhie suggests 

> The highest-quality Anna’s hummingbird genome was achieved by generating haplotype-separated contigs using Pacific Biosciences continuous long reads, followed by scaffolding with 10× Genomics linked reads, Bionano optical maps and Phase Genomics Hi-C; gap filling, base call polishing and manual curation completed the workflow.



So haplotipe resolves contigs assembly > Scaffolding with 10x genomics > Bionano optical map > Hi-C > Gap filling > Manual correction







Watching SMRT Leiden 2021 | Leaf to go presentation By Luca Ermini


[LongQC tool for quality control](https://www.g3journal.org/content/10/4/1193#:~:text=LongQC%20is%20a%20computationally%20efficient,introduced%20at%20the%20sequencing%20stage.)  

Utilise cela pour ton article


[svpack]() | new hifi asssembler called [LJA](https://www.biorxiv.org/content/10.1101/2020.12.10.420448v1.full)



# Make a gene density file using RIdeogram R package


```r

require(RIdeogram)
library(RIdeogram)

Mihong_karyotype <- read.table("Mihong_Chr.txt", sep = "", header = T, stringsAsFactors = F)
Mihong_karyotype


Mihong_gene_density-2 <- GFFex(input = "Mihong_GFF.gff", karyotype = "Mihong_Chr.txt", feature = "gene", window = 100000)
Mihong_gene_density

ideogram(karyotype = Mihong_karyotype, overlaid = Mihong_gene_density)
convertSVG("chromosome.svg", device = "png")

```


Note:

The karyotype should be formatted as follow:

- Save in txt file in excel. csv doesn't work!!!!!!!
- for the GFF file ........










[A pseudomolecule-scale genome assembly of the liverwort Marchantia polymorpha](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.14602)

Tres bon boulot surout pour la detection des centromeres 

Seydina Issa. Diop,



Oliver Subotic,Alejandro Giraldo-Fonseca,Manuel Waller,Alexander Kirbis,Anna Neubauer,Giacomo Potente,Rachel Murray-Watson,Filip Boskovic,Zoe Bont,Zsofia Hock,Adam C. Payton,Daniël Duijsings,Walter Pirovano,Elena Conti,Ueli Grossniklaus,Stuart F. McDaniel,Péter Szövényi


[CuteSV github](https://github.com/tjiangHIT/cuteSV) | [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02107-y)

Diapo 1

![img](https://github.com/Yedomon/bioinformatics101/blob/main/cuteSV.PNG)

Diapo 2

![img](https://github.com/Yedomon/bioinformatics101/blob/main/cuteSV2.PNG)


Diapo 3

![img](https://github.com/Yedomon/bioinformatics101/blob/main/cuteSV3.PNG)


[Genome Informatics Facility](https://github.com/ISUgenomics)



[Performance of Mapping Approaches for Whole-Genome Bisulfite Sequencing Data in Crop Plants](https://www.frontiersin.org/articles/10.3389/fpls.2020.00176/full?report=reader#h3)

![img](https://www.frontiersin.org/files/Articles/504419/fpls-11-00176-HTML/image_m/fpls-11-00176-g002.jpg)


[Code availability](https://github.com/grehl/benchWGBSmap)



[A perl script to visualize intergated genetic map(s) and synteny map(s) at chromosome level](https://github.com/XuepengSun/MapViewer)

```perl
version: v0.5 (03/09/2018)
Useage:
	
	perl MapViewer.pl [options]
	
options:
	[input]
		-syn	synteny file name [if set two files, "map" is not allowed; <=2]
		-map	map file name [multiple files can be specified with comma]
		-gap	gap file, gaps between scaffods in the pseudochromosome
		-help	this help information
		
	[display]
			parameters for display can be adjusted within the script
			
	[file format]
		synteny file: 
			<chr_id>\t<chr_length>\t<chr_start>\t<chr_end>\t<ref_chr_id>\t<ref_chr_length>\t<ref_chr_start>\t<ref_chr_end>
		
		map file:
			<chr_id>\t<chr_length>\t<chr_pos>\t<LG_id>\t<LG_distance>\t<LG_position>
			
		gap file:
			<chr_id>\t<gap start position>\t<gap length>
			

```

![img](https://github.com/Yedomon/bioinformatics101/blob/main/Chr1.png)


- #### [ISMB 2020 Training on long read transcriptome](https://zenodo.org/record/4732317#.YI1yLu3iuUk)


- #### Just found today this south american bioinfomatics society [github](https://github.com/eead-csic-compbio?tab=repositories)



- #### Variant calling


Amaizing tutorial with GATK from ScilifeLab [NGS workflow](https://scilifelab.github.io/courses/ngsintro/1805/labs/NGS_workflow)

An other one found [here](https://github.com/ReiGao/GWSBE/blob/master/2.bwa.sh)

An other one [here](http://costalab.org/wp-content/uploads/2017/05/lecture_3_ngs.pdf)

Great [PPT1](https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/presentations/ohbm2014/imggen/Nho-ImgGen-WGSeqPractical.pdf0)

Great [PPT2](https://scilifelab.github.io/courses/ngsintro/1509/slides/NGS_AJ.pdf)

[Pipeline explained](https://learn.gencore.bio.nyu.edu/variant-calling/pre-processing/)

Amaizing tools for Hi-C visualisation and other genomics data | available on conda

- #### 2021 | [CoolBox: A flexible toolkit for visual analysis of genomics data](https://www.biorxiv.org/content/10.1101/2021.04.15.439923v1)

- #### 2020 | [Cooler: scalable storage for Hi-C data and other genomically labeled arrays](https://academic.oup.com/bioinformatics/article/36/1/311/5530598)


BREAKING NEWS Today 18 april 2021| Heng Li [tweet](https://twitter.com/lh3lh3/status/1383515558306451465)


> Hifiasm v0.15 released with much improved Hi-C mode especially for non-human species. Also for HiFi-only data, it now outputs two partially phased complete assemblies in addition to primary+alternate contigs. 


So it is possible to get the phased diploid assembly. Then perform a scaffolding later.... NB: Hifiasm is not yet a scaffolder...


Great! 


A good idea is to use hifi data + Hifiasm + hi-c data + SALSA or 3DNA to get quite good diploid accurate chromosome scale data. 







Long-read-tools.org: an interactive catalogue of analysis methods for long-read sequencing data [WEBSITE](https://long-read-tools.org/)



[Samtools markdup for duplicate removal or Picard?](https://www.researchgate.net/post/Samtools_markdup_for_duplicate_removal_or_Picard)


question: transpose-multiple-rows-into-a-single-column


Answer [here](https://www.linuxquestions.org/questions/linux-newbie-8/transpose-multiple-rows-into-a-single-column-873623/)

For space tab delimitade
```pyton 

awk 'RT' RS=" " test.input

```
For comma separeted file

```python

awk 'RT' RS="," test.comma.input > out

```


- #### [snp_CAPS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4785002/pdf/66_244.pdf) | [Lee 1](http://www.koreabreedjournal.org/journal/view.html?doi=10.9787/KJBS.2014.46.2.116) | [Lee 2](http://www.plantbreedbio.org/journal/view.html?volume=8&number=3&spage=293&year=2020)



- #### [Plant Breeding and Genomics Learning Lessons](https://plant-breeding-genomics.extension.org/plant-breeding-and-genomics-learning-lessons/)

- #### [SNP-Based Genetic Linkage Maps for Potato](https://plant-breeding-genomics.extension.org/snp-based-genetic-linkage-maps-for-potato/)

- #### Integrating Molecular Biology and Bioinformatics Education| [Paper](https://www.degruyter.com/document/doi/10.1515/jib-2019-0005/html) | [AppliedGenomeResearch](https://github.com/bpucker/AppliedGenomeResearch)


- #### [Interactive introduction to statistics and probabilities](https://seeing-theory.brown.edu/#firstPage)



- #### [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/)



Blog

- #### [sr-c blog](https://sr-c.github.io/)

- #### [zhugueng blog](https://xuzhougeng.blog.csdn.net/)

- #### [jianshu blog](https://www.jianshu.com/u/740f4b0f11e9)

- #### [ZHUBLOG](https://zhuanlan.zhihu.com/)




Paper


[High-density genetic map using whole-genome resequencing for fine mapping and candidate gene discovery for disease resistance in peanut](https://onlinelibrary.wiley.com/doi/epdf/10.1111/pbi.12930)


Flash


- #### [genome annotation v1](https://blog.csdn.net/u012110870/article/details/82500684) key words `基因组注释植物博客`


- #### [rnaseq nextflow](https://sr-c.github.io/2020/05/25/RNA-seq-DiffExp-analysis/)

- #### [jvciMCscan python](https://sr-c.github.io/2019/01/11/jcvi-MCscan/) | [original tutorial](https://github.com/tanghaibao/jcvi/wiki/MCscan-%28Python-version%29) | [additionnal tuto](https://www.jianshu.com/p/39448b970287)

- #### [MCMC V1](https://www.jianshu.com/p/b12e058c6597) [MCMC V2](http://www.chenlianfu.com/?p=2974)  [MCMC Tree](https://www.jianshu.com/p/b12e058c6597)

- #### [paml parat ka ks codeml](https://blog.csdn.net/weixin_42376118/article/details/112065784) | [ka ks v2](http://blog.sciencenet.cn/blog-3433349-1241328.html)


- #### [pipeline for comparative genomics and evolutionnary](https://blog.csdn.net/qq_36608036/article/details/109466468?utm_medium=distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-2.baidujs&dist_request_id=1328626.667.16153013339855743&depth_1-utm_source=distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-2.baidujs)


- #### [rna SEQ](https://zhuanlan.zhihu.com/p/61847802)

- #### [CAFE tutorial](https://www.jianshu.com/p/146093c91e2b) | [original](https://iu.app.box.com/v/cafetutorial-pdf)
- #### [wgdi](https://xuzhougeng.blog.csdn.net/article/details/114013801)
- #### [HiC pro singularity](https://xuzhougeng.blog.csdn.net/article/details/109892648)

- How to use HPC job and so on [msu edu high coumputer plateform](https://wiki.hpcc.msu.edu/) 
- Bioinfo small [tuto](https://wiki.hpcc.msu.edu/display/ITH/Bioinformatics)
- Example [GTAK4](https://wiki.hpcc.msu.edu/display/ITH/GATK4) 
- install repeatmasker with a modified conda based [strategy](https://xuzhougeng.blog.csdn.net/article/details/102804531)
- [Installing maker with conda](https://wiki.hpcc.msu.edu/display/ITH/Installing+maker+using+conda)
- [CAFE](https://xuzhougeng.blog.csdn.net/article/details/102804514)
- #### ggplot2 [example1](https://www.jianshu.com/p/68aa08ef9f61)
- #### [mercury](https://www.jianshu.com/p/61fefb9a9c5f)
- #### [Phylogenomic_Tutorial || Bayesian Phylogenetic Inference](https://www.jianshu.com/p/7983618c51d4)
- #### [Phylogenomic_Tutorial || ML_Tree inference](https://www.jianshu.com/p/282a94b50418)



