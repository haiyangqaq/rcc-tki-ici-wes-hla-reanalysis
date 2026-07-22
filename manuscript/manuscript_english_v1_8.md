# Treatment-context dependence of mutation burden and HLA class I zygosity in renal cell carcinoma immunotherapy: a multi-cohort reanalysis and allele-level peptide audit

## Abstract

### Background

Mutation burden and human leukocyte antigen class I (HLA-I) diversity describe different components of tumor antigenicity, but their associations with immune checkpoint inhibitor (ICI) outcomes are inconsistent in renal cell carcinoma (RCC). We examined whether a mutation-HLA association observed in a lenvatinib-pembrolizumab cohort was identifiable, reproducible across treated RCC cohorts, and supported by allele-level peptide predictions.

### Methods

We rebuilt a published 24-patient RCC cohort from the official mutation annotation format and HLA supplementary tables and an independently audited progression-free survival (PFS) endpoint. Mutation burden was represented by protein-altering mutation count, not by a recalculated clinical tumor mutational burden. Discovery analyses used Kaplan-Meier estimates, log-rank tests, 12-month restricted mean survival time (RMST), fixed-size permutation, continuous Cox models, and threshold sensitivity analyses. A prespecified mutation-only analysis was performed in 35 nivolumab-treated patients from the Miao/DFCI cohort. Contextual reanalyses included publicly available Braun, JAVELIN Renal 101, Chowell, Nature Cancer, and MSK datasets. The MSK analysis joined ICI outcome data to MSK50K HLA fields for 138 patients and controlled multiplicity across 20 conditional and 10 interaction tests. Finally, all 63 retained position-consistent HLA-A*02:01-versus-partner 9-mer events were rescored with MHCflurry 2.2.1.

### Results

The three HLA-A-homozygous patients in the discovery cohort were all above the official mutation-count median and had PFS events. Their PFS differed from that of the remaining 21 patients (log-rank P=0.00045; fixed-size permutation P=0.0094; 12-month RMST difference, -4.40 months). However, the high-burden/HLA-homozygous composite was identical to HLA-A homozygosity, the low-burden/homozygous cell was empty, and continuous mutation count was not associated with PFS (HR per IQR, 1.48; 95% CI, 0.72-3.05; P=0.283). In Miao/DFCI, higher protein-altering mutation count was associated with shorter PFS (HR per IQR, 2.38; 95% CI, 1.12-5.04; P=0.024), whereas larger nivolumab and avelumab-axitinib cohorts did not reproduce a simple mutation-HLA association. In MSK high-TMB ccRCC, any HLA-I homozygosity was associated with shorter overall survival (HR, 4.83; 95% CI, 1.89-12.33; FDR=0.0196), but the strongest TMB-by-HLA interaction did not survive correction (FDR=0.293). In the complete paired peptide universe, the partner allele won all four prediction metrics in 29/63 events (binomial P=0.615); allocation was instead strongly associated with the C-terminal residue class.

### Conclusions

Mutation burden and HLA-I zygosity showed cohort- and treatment-context-dependent associations in RCC. The discovery composite is not an identifiable mutation-by-HLA interaction, and the MSK result is a post hoc conditional association rather than direct replication. Allele-level predictions support terminal-residue-dependent complementarity, not global rescue by a second HLA-A allele. These findings define a testable boundary for prospective RCC immunogenomic studies but do not establish a clinical biomarker.

**Keywords:** renal cell carcinoma; immune checkpoint inhibitor; tyrosine kinase inhibitor; HLA class I; mutation burden; HLA zygosity; neoantigen prediction

## Introduction

ICI-based combinations are standard first-line treatments for advanced clear cell RCC (ccRCC), including regimens that combine PD-1 blockade with vascular endothelial growth factor receptor-directed tyrosine kinase inhibitors (TKIs) [2,3,22,23]. Durable benefit remains heterogeneous, and biomarkers developed in highly mutated tumors have not transferred reliably to RCC [13]. In particular, total mutation burden has shown little or inconsistent association with outcome in RCC ICI trials [4-8].

Mutation burden and HLA-I genotype represent distinct biological quantities. A somatic mutation may generate a candidate antigenic sequence, whereas HLA genotype constrains which peptides can be presented [20]. Pan-cancer studies have associated maximal HLA-I heterozygosity and greater HLA evolutionary divergence (HED) with ICI outcome [11,12]. However, ccRCC has a modest mutation burden, marked copy-number and microenvironmental heterogeneity, and treatment-specific immune effects [4,18,19]. These features weaken the assumption that more mutations or broader germline HLA diversity should have a uniform clinical direction.

The source study for the present analysis evaluated 24 patients with RCC treated with lenvatinib plus pembrolizumab and reported that mean HED was associated with PFS and duration of response [1]. It also reported four patients homozygous at one or more HLA-I loci and a nonsignificant adverse zygosity trend. Therefore, neither HLA diversity nor the conceptual combination of mutation burden and HLA genotype is novel by itself. The unresolved questions are narrower: whether a locus-specific mutation-HLA pattern is statistically identifiable in this cohort, whether its direction transfers to other treated RCC cohorts, and whether allele-level peptide comparisons support a general presentation rescue mechanism.

We addressed these questions using a source-hierarchical reanalysis. The discovery analysis was rebuilt from the official Lee mutation and HLA tables, with PFS censoring independently checked against the published figures. Supportive and non-supportive external datasets were analyzed or summarized under explicit directness labels. A previously selected 29-event peptide analysis was replaced by a complete 63-event paired universe to remove circular inference. The resulting study is a multi-cohort boundary analysis rather than a validation study.

## Methods

### Study design and evidence hierarchy

This was a secondary analysis of deidentified published data. Evidence was ordered as follows: official study supplements and linked clinical endpoints; independently accessible treated-cohort patient-level data; cohort-level proxy analyses; and biological context without linked outcomes. Discovery, external association, interaction, peptide prediction, and ligandomics results were not treated as interchangeable forms of validation.

### Discovery cohort and source reconstruction

The discovery cohort comprised 24 patients from the published phase Ib/II lenvatinib-pembrolizumab RCC study (NCT02501096; BioProject PRJNA610643) [1]. Somatic variants were obtained from Lee Supplementary Data 1 (`MAF exome sequencing`), and HLA-I genotypes and published mean HED values were obtained from Supplementary Data 2. After deduplication by patient, genomic position, alleles, and variant classification, the official table contained 1,296 protein-altering records. Per-patient protein-altering mutation count was used as the mutation-burden measure. Because the official supplement does not provide a callable megabase denominator, this variable was not labeled TMB.

PFS times and status were taken from the locked discovery survival table. Event coding was reconciled against the published Kaplan-Meier risk table and swimmer plot: 16 patients had progression events and 8 were censored. A historical `Event` field that conflicted with 10 records was excluded. HLA-A, HLA-B, and HLA-C homozygosity were defined by equality of the two published alleles at each locus. Published and internal genotype versions agreed at all three loci for 21/24 patients; all primary discovery analyses used the published genotype.

### Discovery statistical analysis

The primary descriptive state comprised patients with mutation count at or above the cohort median and homozygosity at any HLA-I locus. In this cohort, that state was exactly identical to HLA-A homozygosity. PFS was summarized with Kaplan-Meier estimates and compared by a two-sided log-rank test. A 12-month RMST difference was calculated. A fixed-size permutation diagnostic compared the observed log-rank statistic with all allocations of three patients among 24; because the state was selected after observing the cohort, this permutation does not correct for the full model-selection process.

Mutation count was also modeled continuously per interquartile range (IQR) with a univariable Cox model. Proportional hazards were checked using scaled Schoenfeld residuals. Sensitivity analyses replaced the median threshold with the upper tertile and separately evaluated any-locus homozygosity. No multivariable interaction model was fitted in the discovery cohort because no HLA-A-homozygous patient had a mutation count below the median.

### Miao/DFCI mutation-only analysis

The public cBioPortal `ccrcc_dfci_2019` dataset was linked deterministically at patient level [5]. The analysis included 35 nivolumab-treated ccRCC patients, 27 PFS events, 8 censored observations, and 2,824 protein-altering mutation records. Per-patient protein-altering mutation count was represented as a rank variable and scaled by its IQR to limit sensitivity to count distribution. Cox models were fitted univariably and after adjustment for first-line treatment and log2 nivolumab dose. All three model variables were complete in the fitted 35-patient set; their definitions are given in Supplementary Table S7. Nonlinearity was assessed by adding a quadratic term, and influence was assessed in 35 leave-one-patient-out fits. HLA genotype, HED, HLA loss of heterozygosity (LOH), and patient-specific presentation data were unavailable; this analysis was therefore mutation-only.

### Contextual RCC cohorts

Publicly available RCC datasets were used to define transferability boundaries rather than pooled in a meta-analysis because treatments, endpoints, mutation measures, and HLA constructs differed.

The Braun dataset comprised 261 nivolumab-treated ccRCC patients with linked PFS, overall survival (OS), mutation count, neoantigen count, and HLA-I zygosity [4]. Median-defined high mutation or neoantigen burden and any-locus homozygosity were evaluated using log-rank and Cox analyses. JAVELIN Renal 101 source data provided nonsynonymous variants per megabase, HLA-A subtype count, and PFS for 333 patients treated with avelumab plus axitinib [7]. Within the high-burden stratum, single versus multiple HLA-A subtypes were compared using log-rank and age/TMB-adjusted Cox analyses. The Chowell multi-ICI table included 91 RCC patients with TMB, HED, HLA LOH, PFS, and OS; median-defined high-TMB/low-HED and high-TMB/HLA-LOH groups were evaluated descriptively [14]. Nature Cancer source data provided mean HED by response class for 128 IO/IO-treated patients; Kruskal-Wallis and ER-versus-PD Mann-Whitney tests were used [8]. SNiP-RCC findings were included as published HLA-only context and were not reanalyzed [9].

### MSK RCC ICI analysis

Samstein `tmb_mskcc_2018` ICI outcome records were joined by exact sample identifier to the public `msk_impact_50k_2026` HLA fields [10,15]. One sample per patient was retained, yielding 138 patients with RCC, 51 OS events, and 115 patients classified as ccRCC. TMB was the portal-provided nonsynonymous mutations per megabase. HLA defects included HLA-A homozygosity, any HLA-I homozygosity, HLA-A LOH, any HLA-I LOH, and HED-A=0.

Conditional analyses evaluated each of five HLA definitions within median-high and upper-tertile TMB strata, in all RCC and ccRCC scopes (20 tests). Cox models adjusted for age, sex, and ICI combination versus monotherapy; Benjamini-Hochberg correction was applied across the 20 tests. Ten TMB-by-HLA interaction terms (five continuous and five median-binary) formed a separate correction family. Leave-one-out analyses were performed for the strongest conditional result. MSK50K is distributed under CC BY-NC-ND 4.0; modified patient-level MSK data are therefore not redistributed in the submission package. Reproducibility is provided by study identifiers, code, and aggregate outputs.

### Complete HLA-A*02:01 paired-peptide audit

The archived peptide table contained 63 retained position-consistent, mutation-derived C-terminal 9-mer events from seven HLA-A*02:01-heterozygous patients. This complete retained universe comprised 34 events originally assigned to HLA-A*02:01 and 29 assigned to the partner HLA-A allele. The 29 partner-assigned events were not analyzed alone because doing so conditions inference on the outcome of interest.

Each peptide and its stored N- and C-flanks were rescored against both patient HLA-A alleles using MHCflurry 2.2.1 and presentation model release 2.2.0 [16]. LEN07 was scored with the published HLA-A*34:01 call. A partner win required lower affinity, lower affinity percentile, higher presentation score, and lower presentation percentile. The total partner-win proportion was compared with 0.5 using a two-sided binomial test. C-terminal residues were grouped as basic (K/R/H), aromatic (F/W/Y), aliphatic hydrophobic (L/I/V/M/A), or other. Fisher exact tests compared each class with the remainder, followed by Benjamini-Hochberg correction across four classes. Events were nested within seven patients; these tests are event-level descriptive analyses and are not patient-level inference.

### Ligandomics context

PXD017149 ccRCC HLA-I ligand sequences were summarized as all detected ligands (n=1,457) and tumor-exclusive ligands (n=504) [17]. Position-specific residue frequencies and C-terminal class composition were calculated. These data contain neither discovery-cohort mutations nor linked outcomes and were used only as a real-ligand background.

### Mutational signature quality control

Mutational signatures were not included in the primary model. As a sensitivity audit, PASS biallelic SNVs from Mutect2 filtered VCFs were fitted against COSMIC SBS v3.2 signatures with 100 bootstrap replicates [21]. The three discovery cases highlighted by the historical analysis had zero SBS2 and SBS2+SBS13 detection in all 100 replicates; APOBEC-based labels were therefore removed from the manuscript.

### General statistical considerations

All tests were two-sided. Effect estimates and 95% confidence intervals were prioritized over dichotomous significance. Analyses were conducted in Python 3.13 using pandas, scipy, statsmodels, and lifelines, and in R for the mutational-signature audit. No clinical cutoff was optimized or proposed. The study is exploratory and no sample-size calculation was performed.

## Results

### Official source reconstruction changes the discovery exposure

The official Lee mutation table yielded a median of 48.5 protein-altering mutations per patient (IQR, 27.2), matching the source publication [1]. This measure was not interchangeable with a historical server-derived retained-record count: their patient rankings were weakly correlated (Spearman rho=0.166; P=0.438). The published HLA genotype was therefore paired only with the official mutation count (Figure 1). Three patients were HLA-A homozygous (EIS01, LEN08, and LEN15); one additional patient was homozygous at HLA-B and HLA-C.

![Figure 1](figures/Figure1.png)

### The discovery composite is an early-progression observation but not an identifiable interaction

All three HLA-A-homozygous patients had mutation counts at or above the cohort median, and all had PFS events. Compared with the remaining 21 patients, their PFS distribution differed by log-rank testing (P=0.00045). The fixed-size permutation P value was 0.0094, and the 12-month RMST difference was -4.40 months (Figure 2). These diagnostics show that the three observed PFS records are unusually early relative to other fixed groups of three; they do not account for the post hoc choice of HLA locus, burden threshold, or composite definition.

The composite did not provide evidence that mutation burden modified the HLA effect. It was exactly identical to HLA-A homozygosity, and no patient occupied the low-burden/HLA-A-homozygous cell. Continuous mutation count was not associated with PFS (HR per 27.2-mutation IQR, 1.48; 95% CI, 0.72-3.05; P=0.283). The upper-tertile definition reduced the composite from three patients to two. In contrast, any-locus HLA-I homozygosity alone included four patients and was not associated with PFS (log-rank P=0.253). The data therefore support a locus-specific case observation, not separable mutation and HLA effects.

![Figure 2](figures/Figure2.png)

### Mutation count is associated with shorter PFS in Miao/DFCI but not uniformly across RCC ICI cohorts

In the independent Miao/DFCI cohort, the mutation-rank IQR was 0.486 and the mutation-count IQR was 34 records. Higher mutation rank was associated with shorter PFS (HR per IQR, 2.38; 95% CI, 1.12-5.04; P=0.024), with no evidence of proportional-hazards violation (P=0.834). The raw-count model gave a similar estimate (HR per 34 mutations, 2.43; 95% CI, 1.32-4.48; P=0.004). Adjustment for first-line status and dose retained the association (HR, 2.47; P=0.017). All 35 leave-one-out estimates were greater than 1 (range, 2.04-2.78), and 34/35 had nominal P<0.05. A median split was less informative (log-rank P=0.133), supporting continuous reporting rather than a cutoff (Figure 3).

![Figure 3](figures/Figure3.png)

This direction did not generalize to all external cohorts. In Braun nivolumab-treated ccRCC, neither high mutation burden, high neoantigen burden, HLA-I homozygosity, nor their median-defined combinations were associated with PFS or OS; high mutation burden plus any-locus homozygosity gave P=0.651 for PFS and P=0.608 for OS. In JAVELIN avelumab-axitinib, high-burden patients with one versus multiple HLA-A subtypes had similar PFS (adjusted HR, 0.91; 95% CI, 0.45-1.81; P=0.780). In the Chowell RCC subset, high-TMB/low-HED was not associated with PFS (P=0.332) or OS (P=0.249). Mean HED also did not differ across IO/IO response classes in Nature Cancer source data (Kruskal-Wallis P=0.147; ER versus PD P=0.227). Conversely, SNiP-RCC reported an association between higher HED-B and nivolumab benefit [9].

### The MSK conditional association survives FDR, but the interaction does not

The exact MSK overlap comprised 138 patients and 51 OS events; 115 patients had ccRCC. Among 73 median-high-TMB ccRCC patients, 16 were homozygous at one or more HLA-I loci. Any-locus homozygosity was associated with shorter OS after adjustment for age, sex, and ICI regimen (HR, 4.83; 95% CI, 1.89-12.33; nominal P=0.00098; BH-FDR=0.0196; Figure 4A,B). The proportional-hazards test was nonsignificant (P=0.628), and leave-one-out adjusted HRs ranged from 2.94 to 4.05.

The corresponding effect-modification test was weaker. The binary TMB-high-by-any-HLA-homozygosity interaction had HR 5.44 (95% CI, 1.19-24.90), nominal P=0.029, and BH-FDR=0.293. None of the 10 interaction terms survived correction (Figure 4C). Moreover, the strongest conditional definition was any-locus homozygosity, not the HLA-A-specific definition observed in the Lee cohort, and the endpoint was OS rather than PFS. The MSK finding is therefore a robust post hoc conditional association within that dataset, not a direct replication of the discovery state.

![Figure 4](figures/Figure4.png)

### Complete paired rescoring supports terminal compatibility rather than partner-allele rescue

The complete retained paired universe included 63 events from seven patients. The partner allele satisfied all four superiority criteria in 29/63 events (46.0%), which did not differ from 50% (two-sided binomial P=0.615). Thus, there was no global partner-allele advantage. Recomputed direction agreed with the archived best-allele label in all 63 events, indicating numerical reproducibility of the allocation within the same prediction system.

Allocation was strongly structured by the peptide C terminus (Figure 5). The partner allele won for 24/24 basic-tail events and 3/3 aromatic-tail events, but for only 1/33 aliphatic-hydrophobic events and 1/3 other events. Basic-tail enrichment and aliphatic-hydrophobic depletion remained significant after correction (FDR 1.52×10^-12 and 7.58×10^-14, respectively). These patterns are consistent with allele-specific anchor preferences, but they remain computational because no patient-specific mutant peptide was measured by mass spectrometry or tested in T cells.

![Figure 5](figures/Figure5.png)

The PXD017149 ccRCC ligandome contained basic C-terminal residues in 17.6% of all ligands and 18.5% of tumor-exclusive ligands (Supplementary Figure S5). This confirms that basic termini occur among naturally presented RCC ligands; it does not validate the discovery peptides, patient HLA pairings, or clinical association.

## Discussion

This reanalysis does not confirm the original mutation-HLA bottleneck hypothesis. Instead, it identifies where the hypothesis is and is not supported. In the 24-patient lenvatinib-pembrolizumab cohort, three HLA-A-homozygous patients had relatively high official mutation counts and early PFS events. Yet mutation burden cannot be separated from HLA-A homozygosity because the relevant 2×2 table has an empty cell. A composite label adds no information beyond the three HLA-A-homozygous cases, and its extreme log-rank P value should not be interpreted as validation.

This distinction matters because both parts of the proposed biology have already been studied. Chowell et al. linked HLA-I zygosity and HED to ICI outcome in pan-cancer cohorts [11,12], while Lee et al. reported mean-HED associations in the source cohort and an independent MSK cohort [1]. The present contribution is not rediscovery of HLA diversity. It is the locus-specific reconstruction, explicit demonstration of non-identifiability, cross-cohort falsification, and correction of a circular peptide selection.

The external results argue against a portable RCC mutation-HLA biomarker. Miao/DFCI provided a stable mutation-only association with shorter PFS, and the expanded MSK dataset showed an adverse any-locus homozygosity association within high-TMB ccRCC. However, Braun nivolumab, JAVELIN avelumab-axitinib, Chowell multi-ICI RCC, and Nature Cancer HED analyses were null under related but nonidentical definitions. SNiP-RCC implicated HED-B rather than HLA-A [9]. These differences may reflect treatment, endpoint, disease composition, assay, and variable-definition heterogeneity. They may also indicate that bulk mutation count, clonal neoantigen load, and presentation breadth cannot be substituted for one another. The Nature Cancer association of clonal neoantigens with exceptional IO/IO response, together with B-cell/TLS signals in IO/VEGF response, is consistent with this interpretation [8].

The MSK result deserves particular caution. Its conditional association survived correction across the prespecified 20-test family and was stable to single-patient deletion. Nevertheless, the interaction family did not survive FDR correction. A strong effect within a selected high-TMB stratum does not by itself establish that TMB modifies the effect of HLA zygosity. The MSK cohort also differs from the discovery cohort in treatment and endpoint and uses targeted-panel TMB. It is supportive of further testing, not direct validation.

The complete peptide analysis changes the mechanistic interpretation. Restricting analysis to 29 partner-assigned events necessarily produced apparent partner superiority. When all 63 retained position-consistent events were included, the overall advantage disappeared. The remaining signal was biochemical: basic and aromatic termini favored the partner alleles, whereas aliphatic-hydrophobic termini favored HLA-A*02:01. This is compatible with known allele-specific anchor motifs and shows why scalar labels such as “presentation rescue” can be misleading. It does not demonstrate antigen processing, cell-surface presentation, or T-cell recognition.

Several limitations set the evidence ceiling. The discovery cohort is small and the subgroup was selected post hoc. The original study had already examined HLA diversity, limiting novelty. External datasets use different treatments, endpoints, sequencing assays, and HLA variables, precluding a harmonized pooled effect. The MSK analysis is retrospective and its strongest result is conditional rather than interaction-confirmed. Peptide predictions share a model and event-generation pipeline, are nested within seven patients, and PXD017149 is unlinked biological context. RNA expression, clonality, immunopeptidomics from the study patients, spatial immune features, and functional T-cell assays were unavailable. Finally, the analysis is prognostic within treated cohorts and does not establish treatment-predictive utility without non-ICI comparators.

The appropriate next study is prospective and pre-specified. It should enroll an adequately sized RCC ICI-TKI cohort, retain continuous mutation and clonality measures, genotype all HLA-I loci, predefine zygosity/HED constructs, include an interaction model, and collect RNA, spatial, and peptide-level readouts. Until then, neither the Lee case state nor the MSK conditional result should guide treatment.

## Conclusions

Across RCC immunotherapy datasets, mutation burden and HLA-I zygosity did not define a consistent clinical biomarker. The three-patient Lee observation is real at the record level but statistically inseparable from HLA-A homozygosity, and the MSK high-TMB association lacks multiplicity-corrected interaction support. Complete paired peptide rescoring supports terminal-residue-dependent allele complementarity rather than global partner-allele rescue. The study provides a transparent, testable framework and a map of replication boundaries, not a clinically validated classifier.

## Data availability

Discovery sequencing data are available through BioProject PRJNA610643, with mutation and HLA tables in the Lee et al. supplementary material. The Miao/DFCI cohort is available through cBioPortal study `ccrcc_dfci_2019`. MSK outcome and HLA data are available through cBioPortal studies `tmb_mskcc_2018` and `msk_impact_50k_2026`; use is subject to the source licenses. PXD017149 is available through ProteomeXchange/PRIDE. Aggregate source data and analysis code are included with the submission package. Modified patient-level MSK50K data are not redistributed.

## Code availability

Analysis scripts, a run-order manifest, and aggregate source data accompany this version. Before formal submission, the authors must deposit the frozen code and permitted source data in a permanent public repository and add its DOI here. Modified patient-level MSK records must not be deposited.

## Ethics statement

This study used deidentified public or published secondary data. The corresponding author must confirm the local institutional determination regarding exemption or non-human-subjects status before submission.

## Author contributions, funding, competing interests, and AI disclosure

These declarations require author confirmation and must not be inferred. Complete the author-input template before submission. Computational and language-assistance tools should be disclosed according to the target journal policy; authors remain responsible for all analyses and text.

## Figure legends

**Figure 1 | Source-hierarchical reconstruction of the discovery cohort.** A, Audit design. B, Official protein-altering mutation counts ordered from high to low; orange denotes HLA-A homozygosity and the dashed line denotes the cohort median. C, Locus-specific HLA-I zygosity. D, Patient-level PFS; circles indicate progression events and arrowheads indicate censoring.

**Figure 2 | The discovery state is associated with early PFS but is not an identifiable mutation-HLA interaction.** A, Mutation count by HLA-A zygosity. B, Kaplan-Meier PFS and number at risk for the three HLA-A-homozygous patients versus the remainder. The P value is exploratory and unadjusted for subgroup selection. C, 12-month RMST sensitivity across related definitions. D, Burden-by-zygosity table showing the empty low-burden/homozygous cell.

**Figure 3 | Mutation-only association in the Miao/DFCI nivolumab cohort.** A, Protein-altering mutation counts and PFS status for 35 patients. B, Cox estimates per mutation-rank or count IQR. C, Leave-one-patient-out estimates. D, Median-split Kaplan-Meier curve shown descriptively.

**Figure 4 | Conditional MSK association and cross-cohort heterogeneity.** A, OS in median-high-TMB MSK ccRCC stratified by any-locus HLA-I homozygosity. B, Five adjusted conditional HLA tests within high TMB; FDR values belong to the 20-test family. C, Binary TMB-by-HLA interaction estimates; none survived correction across 10 interaction tests. D, Direction and directness of RCC cohort evidence; rows are not pooled because definitions and endpoints differ.

**Figure 5 | Complete HLA-A*02:01 paired-peptide audit.** A, Partner-minus-HLA-A*02:01 presentation margins for all 63 retained position-consistent events. B, Paired presentation scores by C-terminal class. C, Fraction of partner wins by C-terminal class with four-test BH-FDR. D, Event and partner-win counts by patient.

**Supplementary Figure S1 | Official mutation landscape.** Twenty recurrent or canonical RCC genes from the published Lee mutation table.

**Supplementary Figure S2 | HLA version provenance.** A, Published versus internal genotype agreement. B, Difference between historical and published mean HED values.

**Supplementary Figure S3 | Mutational-signature quality control.** A, PASS biallelic SNV counts. B, SBS2 and SBS2+SBS13 detection rates across 100 strict bootstrap fits.

**Supplementary Figure S4 | Paired prediction metrics.** HLA-A*02:01 versus partner affinity, affinity percentile, presentation score, and presentation percentile for all 63 events.

**Supplementary Figure S5 | PXD017149 ccRCC HLA-ligandome context.** Position-specific C-terminal residue frequencies and terminal-class composition for all and tumor-exclusive ligands. No discovery patient or mutant peptide is represented.

**Supplementary Figure S6 | Full MSK sensitivity family.** Twenty conditional tests and ten interaction tests used for multiplicity control.

## References

1. Lee CH, DiNatale RG, Chowell D, et al. High response rate and durability driven by HLA genetic diversity in patients with kidney cancer treated with lenvatinib and pembrolizumab. Mol Cancer Res. 2021;19:1510-1521. doi:10.1158/1541-7786.MCR-21-0053
2. Motzer R, Alekseev B, Rha SY, et al. Lenvatinib plus pembrolizumab or everolimus for advanced renal cell carcinoma. N Engl J Med. 2021;384:1289-1300. doi:10.1056/NEJMoa2035716
3. Powles T, Albiges L, Bex A, et al. Renal cell carcinoma: ESMO Clinical Practice Guideline for diagnosis, treatment and follow-up. Ann Oncol. 2024;35:692-706. doi:10.1016/j.annonc.2024.05.537
4. Braun DA, Hou Y, Bakouny Z, et al. Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma. Nat Med. 2020;26:909-918. doi:10.1038/s41591-020-0839-y
5. Miao D, Margolis CA, Gao W, et al. Genomic correlates of response to immune checkpoint therapies in clear cell renal cell carcinoma. Science. 2018;359:801-806. doi:10.1126/science.aan5951
6. McDermott DF, Huseni MA, Atkins MB, et al. Clinical activity and molecular correlates of response to atezolizumab alone or in combination with bevacizumab versus sunitinib in renal cell carcinoma. Nat Med. 2018;24:749-757. doi:10.1038/s41591-018-0053-3
7. Choueiri TK, Donahue AC, Braun DA, et al. Integrative analyses of tumor and peripheral biomarkers in the treatment of advanced renal cell carcinoma. Cancer Discov. 2024;14:406-423. doi:10.1158/2159-8290.CD-23-0680
8. Jammihal T, Saliby RM, Labaki C, et al. Immunogenomic determinants of exceptional response to immune checkpoint inhibition in renal cell carcinoma. Nat Cancer. 2025;6:372-384. doi:10.1038/s43018-024-00896-w
9. Tanegashima T, Shiota M, Fujiyama N, et al. Effect of HLA genotype on anti-PD-1 antibody treatment for advanced renal cell carcinoma in the SNiP-RCC study. J Immunol. 2024;213:23-28. doi:10.4049/jimmunol.2300308
10. Samstein RM, Lee CH, Shoushtari AN, et al. Tumor mutational load predicts survival after immunotherapy across multiple cancer types. Nat Genet. 2019;51:202-206. doi:10.1038/s41588-018-0312-8
11. Chowell D, Morris LGT, Grigg CM, et al. Patient HLA class I genotype influences cancer response to checkpoint blockade immunotherapy. Science. 2018;359:582-587. doi:10.1126/science.aao4572
12. Chowell D, Krishna C, Pierini F, et al. Evolutionary divergence of HLA class I genotype impacts efficacy of cancer immunotherapy. Nat Med. 2019;25:1715-1720. doi:10.1038/s41591-019-0639-4
13. Litchfield K, Reading JL, Puttick C, et al. Meta-analysis of tumor- and T cell-intrinsic mechanisms of sensitization to checkpoint inhibition. Cell. 2021;184:596-614.e14. doi:10.1016/j.cell.2021.01.002
14. Chowell D, Yoo SK, Valero C, et al. Improved prediction of immune checkpoint blockade efficacy across multiple cancer types. Nat Biotechnol. 2022;40:499-506. doi:10.1038/s41587-021-01070-8
15. Bandlamudi C, et al. Cancer type-specific variation in patterns of driver alterations across 50,000 tumors. Cancer Cell. 2026. doi:10.1016/j.ccell.2026.03.003
16. O'Donnell TJ, Rubinsteyn A, Laserson U. MHCflurry 2.0: improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. Cell Syst. 2020;11:42-48.e7. doi:10.1016/j.cels.2020.06.010
17. Reustle A, et al. Integrative omics and HLA-ligandomics analysis to identify novel drug targets for ccRCC immunotherapy. Genome Med. 2020;12:32. doi:10.1186/s13073-020-00731-8
18. Kinget L, Naulaerts S, Govaerts J, et al. A spatial architecture-embedding HLA signature to predict clinical response to immunotherapy in renal cell carcinoma. Nat Med. 2024;30:1667-1679. doi:10.1038/s41591-024-02978-9
19. Cancer Genome Atlas Research Network. Comprehensive molecular characterization of clear cell renal cell carcinoma. Nature. 2013;499:43-49. doi:10.1038/nature12222
20. Schumacher TN, Schreiber RD. Neoantigens in cancer immunotherapy. Science. 2015;348:69-74. doi:10.1126/science.aaa4971
21. Alexandrov LB, Kim J, Haradhvala NJ, et al. The repertoire of mutational signatures in human cancer. Nature. 2020;578:94-101. doi:10.1038/s41586-020-1943-3
22. Rini BI, Plimack ER, Stus V, et al. Pembrolizumab plus axitinib versus sunitinib for advanced renal-cell carcinoma. N Engl J Med. 2019;380:1116-1127. doi:10.1056/NEJMoa1816714
23. Motzer RJ, Tannir NM, McDermott DF, et al. Nivolumab plus ipilimumab versus sunitinib in advanced renal-cell carcinoma. N Engl J Med. 2018;378:1277-1290. doi:10.1056/NEJMoa1712126
