# bioinformatics101

- #### Variant calling


Amaizing tutorial with GATK from ScilifeLab [NGS workflow](https://scilifelab.github.io/courses/ngsintro/1805/labs/NGS_workflow)




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



