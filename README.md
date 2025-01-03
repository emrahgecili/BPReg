# BPReg
Bayesian penalized regressions for non-stationary Gaussian linear mixed effects model.

Gecili, E, Sivaganesan, S, Asar, O, Clancy, JP, Ziady, A, Szczesniak, RD. Bayesian regularization for a nonstationary Gaussian linear mixed effects model. 
Statistics in Medicine. 2022; 41( 4): 681– 697. doi:10.1002/sim.9279

Description: Bayesian penalized regressions for non-stationary Gaussian linear mixed effects model (Diggle et al. 2015). 
Specifically, it assumes integrated brownian motion for W as stochastic process. Returns posterior estimates and credible intervals.
It additionally contains functions for diagnostic checks of posteriro samples such as tracing plots,etc.



# Bayesian regularization for a nonstationary Gaussian linear mixed effects model

In omics experiments, estimation and variable selection can involve thousands of proteins/genes observed from a relatively small number of subjects. Many regression regularization procedures have been developed for estimation and variable selection in such high-dimensional problems. However, approaches have predominantly focused on linear regression models that ignore correlation arising from long sequences of repeated measurements on the outcome. Our work is motivated by the need to identify proteomic biomarkers that improve the prediction of rapid lung-function decline for individuals with cystic fibrosis (CF) lung disease. We extend four Bayesian penalized regression approaches for a Gaussian linear mixed effects model with nonstationary covariance structure to account for the complicated structure of longitudinal lung function data while simultaneously estimating unknown parameters and selecting important protein isoforms to improve predictive performance. Different types of shrinkage priors are evaluated to induce variable selection in a fully Bayesian framework. The approaches are studied with simulations. We apply the proposed method to real proteomics and lung-function outcome data from our motivating CF study, identifying a set of relevant clinical/demographic predictors and a proteomic biomarker for rapid decline of lung function. We also illustrate the methods on CD4 yeast cell-cycle genomic data, confirming that the proposed method identifies transcription factors that have been highlighted in the literature for their importance as cell cycle transcription factors.

# Method

Our methodology and details of prediction model are described in our main article [1]. The
original article provides additional real data applications with a real cystic-fibrosis
(CF) lung function and proteomic datasets and simulations for different data settings.

# CD4 Yeast Cell-Cycle Genomic Data Analysis

We analyzed a subset of the yeast cell-cycle gene expression, which has been previously considered by others [2-4].These data were longitudinally collected in the CDC15 experiment performed by Spellman et al. [5] where genome-wide mRNA levels of 6178 yeast open reading frames in a two cell-cycle period were measured at M/G1-G1-S-G2-M stages. It is critical to identify transcription factors (TFs) that regulate the gene expression levels of cell cycle-regulated genes to better understand the mechanism underlying cell-cycle process. The subset data that we analyze in this application are obtained from PGEE R package [6], and it consists of 297 cell-cycle regularized genes observed over 4 time points at G1 stage and the standardized binding probabilities of a total of 
 TFs obtained from a mixture model approach of Wang, Chen and Li [4]. The response variable is the log-transformed gene expression levels and the covariates are the binding probabilities of 96 TFs. The proposed penalized regression models were applied to this data to identify the TFs that influence the gene expression level at G1 stage of the cell process. Again only $\beta_{k}$'s, the coefficients of TFs, are penalized to induce variable selection.

In our analysis, BL identified 3 TPs (MBP1, STB1, and NDD1), BAL identified 4 TPs (MBP1, STB1, NDD1, and MTH1), BEN identified 5 TFs (MBP1, STB1, NDD1, SWI6, and FKH2), while BR were able to picked 14 TPs (ABF1, FKH2, GRF10.Pho2., HIR1, HIR2, MET4, MTH1, NDD1, SWI4, YAP6, CIN5, HSF1, MCM1, and SKN7). The computation time for BEN was 6 hours 45 minutes for this specific example when MCMC size was 20k which was long enough for convergence. Although the subjects in this study do not have long sequences of repeated measurements, the proposed methods were able to identify important TFs that have already been verified by some biological experiments using genome-wide binding techniques. For example, MBP1 is a crucial transcription factor involved in cell cycle progression from G1 to S stage; NDD1 regulate G2/M genes through binding to their promotes; function of FKH2 is activation of its M stage-specific target genes and it is a cell cycle activator for genes in GFKH2 during the G2 stage; STB1 encodes a protein that contributes to the regulation of SBF and MBF target genes; expression is cell cycle-regulated in stages G1 and S. BR resulted in more discoveries and all of these additional TFs have been reported as key cell cycle TFs in different stages. We refer to studies by [2-4] for additional details on context of the TFs.

We incorporated a variety of Bayesian shrinkage approaches to perform variable selection for a Gaussian linear mixed effects model with nonstationary covariance structure, thereby incorporating the complicated structure of long sequences of repeated, mistimed measurements on the response variable. For the genomic data example, our approaches corroborated TFs that have been identified previously in the literature.

![Posterior mean estimates and their 95\% credible intervals for the coefficients for select TFs from our four models.](https://github.com/emrahgecili/BPReg/blob/master/forest%20plots%20for%20all%20four%20models.png)


# Example R code for BEN model

library(mcmcplots)

library(BPReg)

library(PGEE)

data(yeastG1) 

data <- yeastG1 

pros=data[,4:99]

x=as.matrix(pros)

y=data$y

id=data$id 

time<-data$time

cova<-cbind(x)

x<-scale(cova)

MCMCsamp=BPReg(y,x,time,id,n.samples =8000,r.init = 1,delta.init = 8,method = "BEN")

M=MCMC_aburn(burnin=7000)

posterior_summary(M, probs = c(0.025, 0.975),robust =FALSE)

colnames(M)[2]<-"visit_age"

MCMCplot(M[,c(2:28)],params=colnames(M)[2:28],ref_ovl = TRUE,main="Parameter estimates with BEN",xlab=NULL)

Inf_criteria(M)

MCMCtrace(M, params =colnames(M),ISB = FALSE,pdf = FALSE)
mcmcplot(M)

# References

[1] Gecili E, Sivaganesan S, Asar O, Clancy JP, Ziady A, Szczesniak RD. Bayesian regularization 
for a nonstationary Gaussian linear mixed effects model. Stat Med. 2022 Feb 20;41(4):681-697. 
doi: 10.1002/sim.9279. Epub 2021 Dec 12. PMID: 34897771; PMCID: PMC8795479.

[2] Banerjee N, Zhang MQ. Identifying cooperativity among transcription factors controlling
the cell cycle in yeast. Nucleic Acids Res. 2003;31(23):7024-7031. doi:10.1093/nar/gkg894

[3] Tsai HK, Lu HH, Li WH. Statistical methods for identifying yeast cell cycle transcription
factors. Proc Natl Acad Sci U S A. 2005;102(38):13532-13537. doi:10.1073/pnas.0505874102

[4] Wang L, Chen G, Li H. Group SCAD regression analysis for microarray time course
gene expression data. Bioinformatics. 2007;23(12):1486-1494. doi:10.1093/bioinformatics/btm125

[5] Spellman PT, Sherlock G, Zhang MQ, et al. Comprehensive identification of cell cycle-regulated 
genes of the yeast saccharomyces cerevisiae by microarray hybridization. Mol Biol Cell. 1998; 9(12): 3273-3297.

[6] Inan G, Wang L. PGEE: an R package for analysis of longitudinal data with high-dimensional covariates.
R J. 2017; 1(9): 393-402. doi:10.32614/RJ-2017-030
