# A multimodal simulation framework for benchmarking lineage recording, clonal reconstruction, and cell-state inference

Aidan Cook^1,2^, Aaron McKenna^1,2,3^

^1^ Graduate Program in Quantitative Biomedical Sciences, Dartmouth College, Hanover, NH  
^2^ Department of Molecular and Systems Biology, Dartmouth College, Hanover, NH  
^3^ Dartmouth Cancer Center, Dartmouth College, Lebanon, NH

## Abstract

Recent advances in lineage tracing have enabled increasingly high-resolution reconstruction of developmental and tumor phylogenies. Prospective CRISPR recorders and endogenous mitochondrial mutations capture complementary aspects of cellular history, but their relative and joint performance depends on recorder design, mutation rate, molecular recovery, mitochondrial drift, population dynamics, and downstream inference. We developed a modular simulation framework that jointly models evolving mitochondrial genomes, engineered lineage barcodes, cell-type transitions, cell division, and selective cell loss. The framework generates ground-truth trees and observation matrices for benchmarking lineage reconstruction and clade-level cell-state inference across experimental parameter spaces. In addition to its native forward-time population simulator, the software can import lineage and division outputs from the three-dimensional agent-based simulator PhysiCell, convert retained-parent identifiers into an event-resolved binary tree, and replay base-editing and mitochondrial recording along measured branch durations. A reference-conditioned scDesign3 layer then generates transcriptomes for the same terminal cells using cell type, lineage depth, spatial position, and optional recording-derived covariates. Together, these components provide a bridge from abstract recorder design to spatial tumor growth and multimodal single-cell readout. We use this framework to evaluate when engineered and endogenous lineage signals recover global or subclonal topology and to identify assumptions that limit interpretation of lineage-linked cell-state analyses.

## Introduction

The mitochondrial genome has drawn significant interest for lineage tracing. Cells contain hundreds to thousands of copies of the 16.6-kb circular mitochondrial genome, which accumulates mutations at substantially higher rates than nuclear DNA. As these mutations arise and are inherited during cell division, they generate heteroplasmy, the coexistence of multiple mitochondrial alleles within a single cell. Although heteroplasmy has long been studied for its effects on cellular fitness and disease, these naturally occurring variants also provide a rich source of divergent markers for reconstructing cellular lineage histories. Mitochondria also undergo random intracellular drift, fusion and fission events, non-Mendelian inheritance patterns, and mitophagy (Elson et al. 2001; Stewart and Chinnery 2015; Wallace and Chalkia 2013; W. Chen, Zhao, and Li 2023; Tilokani et al. 2018). These factors, coupled with non-Mendelian transmission patterns of mitochondria from parent to daughter cells, raise questions about the utility of mitochondrial mutational signatures in reconstructing lineage relationships (Li et al. 2026).
Meanwhile, prospective lineage tracing methods have improved the power of lineage reconstruction by directly marking cells at a given timepoint and tracking their descendants longitudinally (McKenna and Gagnon 2019). These approaches are sufficiently powered to study cell fate, the effects of intrinsic and extrinsic perturbations on development, and clonal and subclonal relationships between cells (Wagner and Klein 2020). Indeed, this enhanced resolution comes at the expense of a high experimental burden to successfully engineer the diverse components of these systems, adding to the demand for developing tools that ensure experiments are well designed before beginning them.

Among the most popular and powerful prospective lineage tracing approaches are CRISPR-Cas systems in which Cas/gRNA complexes make successive edits in synthetic nucleotide barcodes that have been engineered into an in vitro or in vivo system. Early iterations of this technology relied on a Cas9 nuclease to generate double-stranded breaks (DSBs) at predefined genomic loci; error-prone repair mechanisms then generate heritable insertions and deletions at these target sites due to the stochastic outcomes of processes such as non-homologous end joining (NHEJ) (McKenna et al. 2016; Spanjaard et al. 2018; Frieda et al. 2017). Mutations resulting from these repair processes serve as heritable markers that ultimately define sublineages. A key limitation of nuclease-based approaches is the dropout of intervening barcode sequences when simultaneous deletion events occur at nearby target sites, which can yield decreased edit diversity at a given target site across cells (Salvador-Martinez et al. 2019).

To address these limitations, more recent iterations of CRISPR-Cas-based lineage tracing approaches have been developed by modifying engineered components of the lineage tracing system (Anzalone, Koblan, and Liu 2020). One particularly promising emerging technology is CRISPR-Cas base editing systems. These systems rely on the fusion of a catalytically inactive endonuclease and a deaminase, circumventing the undesired byproducts of NHEJ by inducing predictable point mutations at nucleotide targets within a specified editing window (Porto et al. 2020; Komor et al. 2016). Notably, since DSBs are not generated in the base conversion process, more predictable mutational patterns are possible. Nuclease systems involving homing guide RNAs (hgRNAs) allow for iterative edits of their own spacer sequences by incorporating protospacer adjacent motif (PAM) sites within each guide locus, thus allowing for the evolution of diverse nuclease outcomes through self-targeting (Kalhor et al. 2018). Prime editing via prime editing guide RNAs (pegRNAs) offers yet another alternative to generate heterogeneous mutation profiles without double stranded breaks (Anzalone et al. 2019; P. J. Chen and Liu 2023). The key component of this approach is a fusion protein between a programmable nickase, an engineered reverse transcriptase, and a pegRNA that contains a sequence corresponding to the desired edit to be made at a programmed target site.

Previous in silico works have explored simulation frameworks that independently model the reconstructive power of specific CRISPR-Cas lineage tracing systems or retrospective modalities (Liu, Zhang, and Yang 2024; Jones et al. 2020). However, experimental design also depends on how recording modalities interact with population structure, spatial growth, cell state, and molecular sampling. A useful simulator should therefore separate the process that generates a cell lineage from the processes that write and recover molecular lineage information.

Here we introduce a framework with two complementary lineage-generation modes. Its native forward-time simulator jointly models population dynamics, cell-state transitions, mitochondrial inheritance, and engineered barcode editing. Its fixed-tree mode instead accepts division times and topology from an external simulator and overlays molecular recording on those branches. We demonstrate the latter with PhysiCell, an agent-based framework for three-dimensional multicellular systems (Ghaffarizadeh et al. 2018). We further connect terminal cells to reference-conditioned transcriptomes generated with scDesign3 (Song et al. 2024), producing aligned lineage, spatial, recording, and expression outputs. We use these capabilities to explore when prospective and retrospective signals recover accurate single-cell phylogenies and when they support downstream clonal fate inference.

## Results

### Framework overview: native population simulation and fixed-tree replay

The framework separates three layers that are often coupled in lineage-tracing experiments: lineage generation, molecular recording, and terminal-cell observation. In native mode, all three layers are simulated together. A founder population expands according to stochastic cell-cycle schedules, optional cell-type transitions, and cell loss; mitochondrial and engineered-barcode profiles evolve through the same divisions. In fixed-tree mode, lineage topology and division times are supplied externally, while molecular profiles are generated by replaying the recording model along the imported branch durations. Both modes emit ground-truth trees and terminal-cell mutation profiles that can be processed through the same reconstruction and clade-analysis workflows.

For spatial tumor simulations, a lineage-enabled PhysiCell project grows a three-dimensional tumor from a founder cell to a user-defined population threshold. The exported persistent-parent division log is normalized into an event-resolved binary tree, after which base-editing and mitochondrial recording are simulated for every branch. Terminal PhysiCell identifiers are retained as shared sample keys across Newick trees, barcode matrices, mitochondrial heteroplasmy tables, and optional transcriptomes. This architecture enables recorder designs to be compared on either abstract cell populations or externally generated, spatially explicit lineages without changing the downstream observation model.

### Prospective lineage-recorder performance

We first benchmark the framework by simulating two prospective lineage-tracing systems: a Cas9 nuclease and an adenine base editor, parameterized by observed edit counts in FLARE and BASELINE, respectively. Three replicates of each experiment were run, and accuracy scores were mean-aggregated.

Generally, the base editing system outperformed the nuclease system according to tree reconstruction accuracy, here defined as 1 minus normalized Robinson-Foulds distance between prospective barcode-estimated trees and ground truth lineage structure. At initial founder integrations above 20 and perfect integration recovery probabilities, BASELINE was able to perfectly reconstruct underlying tree structure across all cell sampling fractions tested once the cell population had reached sufficient size. While also improving in accuracy with increasing cell sampling fraction and number of founder integrations, FLARE was not able to reliably capture this perfect tree structure.

> **Analysis note:** Resolve the non-monotonic accuracy observed for some BASELINE population-size comparisons before finalizing this subsection.

### Mitochondrial lineage signal and clade recovery

We next use the framework to compare the relative information capacity of mitochondrial lineage tracing across mutation, heteroplasmy, fusion/fission, and inheritance parameterizations. Founder mitochondrial genomes contain low-frequency variants whose allele fractions are drawn from a beta distribution. Terminal mutation profiles are converted to allelic-fraction or thresholded binary matrices and evaluated either as phylogenetic characters or as mitochondrially defined mutation-profile clades.

The progressive bottlenecking analyses yielded the strongest coarse-clade overlap in the current draft results, whereas less constrained inheritance regimes produced more callable variants but weaker agreement with ground-truth clades. These results suggest that marker abundance alone does not determine reconstructive value: transmission and drift determine whether variants remain informative at the scale of interest.

> **Implementation/reanalysis note:** The maintained native simulator currently implements random mitochondrial inheritance. Results labeled as biased, bottleneck, or clustered Models J–L were generated during exploratory model development and should be rerun against a restored, tested implementation or moved to a clearly labeled conceptual analysis before submission. The PhysiCell replay separately implements a fixed-size genome bottleneck, but it does not reproduce organelle-level fusion/fission or variable copy number.

### Clade-level inference of cell-state transition regimes

We next asked whether lineage-defined clades could distinguish two cell-state transition regimes that produce similar global cell-type proportions. Ground-truth topological clades and mitochondrial mutation-profile clades were passed to an expectation-maximization model that estimates the probability that each clade arose under an induced rather than uninduced transition matrix. In the draft analyses, mitochondrial clades performed below topology-defined gold-standard clades, although the bottlenecked model retained the strongest relative signal. Sensitivity analysis suggested that the classifier required the correct direction of the transition bias but was comparatively robust to uncertainty in its exact magnitude.

> **Analysis note:** The 20-clade ground-truth analysis contained a known bug and should not be reported until rerun. Retain five-clade results as provisional.

### Spatial tumor lineage replay with PhysiCell

The fixed-tree workflow was validated using the lineage-enabled `tumor_3D_lineage` PhysiCell project. The end-to-end driver stages the project without modifying the PhysiCell checkout, grows the tumor to a configurable cell-count threshold, exports the complete division log and final-cell table, and replays both recording modalities. An eight-cell integration test produced seven division events, a 15-node event-resolved tree, and matched barcode and mitochondrial outputs for all eight terminal cells. The production workflow is configured for 10,000 cells; full-scale biological benchmarking remains to be performed.

This integration makes it possible to ask whether recorder performance changes with spatially generated branch times, neighborhood structure, and tumor geometry. The current PhysiCell export includes terminal positions and neighbor identifiers but does not yet include per-cell phenotype, substrate exposure, cell-cycle state, or death time. These fields define natural extensions for future spatially conditioned analyses.

### Covariate-linked transcriptome generation

To connect simulated lineage history to a single-cell molecular phenotype, we fit scDesign3 marginal and gene-correlation models to a user-provided reference dataset and generate expression counts for the terminal PhysiCell cells. The target covariate table contains cell type, lineage depth, normalized lineage pseudotime, branch timing, x/y/z position, tumor radius, neighbor count, barcode edit burden, and mitochondrial heteroplasmy summaries. Only covariates represented in the reference and selected in the scDesign3 formula influence expression; all fields remain available as terminal-cell metadata.

The implementation was validated by fitting a reference-conditioned model with lineage pseudotime and barcode edit fraction and by separately fitting a spatial model. Synthetic expression columns exactly matched the terminal sample identifiers used in the recording and tree outputs. This standard integration generates conditionally independent cells given their covariates; it does not impose additional transcriptomic covariance between close relatives.

## Methods

### Software architecture and simulation modes

The software exposes a native forward-time simulator and a fixed-tree recording simulator. The native mode generates population history, molecular profiles, and cell states jointly. The fixed-tree mode imports lineage topology and branch durations from PhysiCell and propagates recording profiles along that supplied history. A shared post-processing layer converts terminal profiles into sequence and character matrices, and an optional scDesign3 stage generates aligned expression counts.

### Native generation of the cell population
Here we describe the processes underlying a simulation in which we are interested in tree reconstruction from both retrospective and prospective mutations.

The native simulator initializes a cell population C(0) = {c0} at time t = 0. Each cell c in C(t) is characterized by the state vector c = (M_mt(c), M_bc(c), τ(c), θ(c)), where:

- M_mt(c) is the sparse mitochondrial mutation matrix of a cell c, a g-by-k matrix, with g mitochondrial genome copies of length k each,
- M_bc(c) is the sparse barcode mutation matrix of a cell c, an m-by-ℓ matrix, with a maximum of m barcode integrations of length ℓ each,
- τ(c) denotes the cell type of cell c, and
- θ(c) denotes the cell- and cell type-specific parameters of cell c.

Let n_mito denote the number of mitochondria in the founder cell. Define a map φ that assigns each genome row i of the g-by-k mutation matrix M_mt(c) to a mitochondrion φ(i), where φ maps {1,...,g} to {1,...,n_mito}. Let the inverse mapping from mitochondria to genome row indices be represented by φ^-1.

Let K be the number of distinct cell types. The discrete cell type τ(c) of a cell is defined when the cell arises from a division and may change to any other cell type at its own division point according to the K-by-K transition matrix Phi0. Binary inducible cell differentiation conditions can be captured by providing an induction timepoint t_I such that all cell type transitions at timepoints t ≥ t_I use the K-by-K transition matrix Phi1, where each row of Phi_r sums to 1 (r = 0 or 1). For a cell c born at t_(c,birth) with type τ(c) = j, its cell cycle length is drawn from an exponential distribution with rate λ_k = 1/μ_k, where μ_k is the mean cell cycle length for cell type k. The next eligible division timepoint for this cell, t_(c,elig), is defined as t_(c,birth) + (cell cycle length). The simulation progresses from t = 0 to t_sim with discrete time points separated by a time increment Δt.

At each t, up to three processes can occur for each cell. For those cells past their eligible division timepoint, cell division occurs; stochastic and heteroplasmy-driven death events prune the population; and all cells alive at t undergo mutation processes that update M_mt and M_bc. During a cell division, daughter cells receive identical copies of parental barcode mutation matrices but perturbed copies of parental mitochondrial mutation matrices. Stochastic death processes occur with probabilities p_death(c) that are defined by τ(c). Transmission of mitochondrial profiles, mutation profile updates, and heteroplasmy-driven death are described in more detail below.

The output of the simulation is a discrete cell type-labeled heterogeneous cell population whose mutation matrices can be further processed to infer lineage relationships. By performing multiple post-hoc processing steps on the same mutation matrices, we can more easily draw comparisons between tree reconstruction performance metrics across a range of possible experimental limitations.

### Mutation-profile representation and evolution

Mutation matrices M_mt and M_bc store numeric representations of site- and integration-/genome-specific deviations from wild-type or baseline sequences across the lineage barcodes of a cell. We assume all barcode target sites are unedited at time t = 0, such that M_bc(c) is the all-zero m-by-l matrix; initial mitochondrial heteroplasmy makes this untrue with mitochondrial barcodes. For efficient representation of mutation profiles, we map nucleotides (A, G, C, T) to integers 1, 2, 3, and 4 respectively.

At each t, the barcode mutational matrices of all living cells undergo uniform and targeted mutational processes in two separate steps; the mitochondrial mutational matrices undergo only the uniform mutational process. Targeted mutational processes differ from their uniform counterparts by iteratively mutating barcode positions that have been predefined as targets and labeled with a discretized High, Medium, or Low edit rate classification to capture the variability in editing efficacy across guide/target combinations. The placement of target sites across the genomic barcode can be randomized or manually specified. Further heterogeneity within edit rate class is introduced by drawing from different quantile-defined regions of the underlying gamma distribution that specifies all target mutation probabilities of a given cell type and mutation type. These mutation types are defined as transitions, transversions, insertions, and deletions. In all mutation processes, we assume that multiple mutations cannot occur at the exact same base position in the same integration or genome, and we allow mutation rates to vary by cell type.

More formally, for each genomic position p and for each mutation type v, we generate the number of new integrations or genomes mutated, E_(p,v), as a binomial draw with N_p(unedited) trials and success probability π_(p,v), where N_p(unedited) is the number of unedited integrations or genomes at position p and π_(p,v) is the probability of an edit of type v occurring in any one genome or integration at p. Note that in uniform editing processes, π_(p,v) = π_v (independent of position). For each v, the E_(p,v) integrations or genomes to be edited are then randomly sampled from the unedited pool at position p.

Inducible editing systems can be modeled in a binary way where prior to induction, all targets and non-targets alike undergo only background mutational processes; after induction, targets selectively undergo an additional round of non-uniform editing at each simulation timepoint.

For each cell at each simulation step, newly acquired mutations are encoded in temporary timepoint-specific mutation matrices R_mt(c_t) and R_bc(c_t). Without loss of generality in the mitochondrial context, a point mutation resulting in nucleotide base B that occurs in integration i at position j is encoded as R_bc(c_t)[i,j] = (integer code for B). Insertions of length L_insert in integration i that begin at position j are encoded as fractional values less than 1, of the form 0.d1 d2 ... dL, where each digit dℓ is the integer code for nucleotide B at that position of the insertion. Deletions of length d that start at position j-d+1 in integration i are represented by the value -1 in each of the d consecutive entries from position j-d+1 through position j.

Let M_bc(c_t) be the incoming barcode mutation matrix, and let M_bc(c_(t+Δt)) be the barcode mutation matrix following these mutational processes, generated by M_bc(c_(t+Δt)) = M_bc(c_t) + R_bc(c_t) (similarly, M_mt(c_(t+Δt)) = M_mt(c_t) + R_mt(c_t)). This algebraic step ensures that mutations accumulated at time t are efficiently encoded into mutational profiles that can again be modified at future time steps.

### Prospective CRISPR recorder models

Each barcode integration is assumed to have an identifiable tag that enables integrations to be compared across cells that descend from the same founder. For each specified target in a CRISPR-Cas system, editing windows may be built from optionally decaying edit rates with increasing distance from target sites. By indicating the size and behavior of this window, the user can allow for or prevent multiple edits occurring within the same target region.

We identify combinations of parameters to emulate the biological mechanisms and editing patterns of CRISPR/Cas prospective lineage tracing systems. A shorthand form of the key parameterizations of each system is given in Table [XXXX].

Conventional Cas9 systems are encoded with a decaying editing window surrounding each target, as well as elevated insertion and deletion rates. Upon mutation of a target site within this editing window, further modifications are prohibited. To mimic the observed behavior of inter-target dropout that results from simultaneous deletions occurring at nearby genomic regions, we encode a dropout radius that probabilistically deletes all intervening genomic sequences between two sites when simultaneous deletions occur.

Base editing systems are modeled by increasing the frequency of either transitions or transversions at target sites, depending on the specified nucleotide conversion pattern. Target windows in a base editing context identify off-target nucleotides with the same identity as the target and encode them as possible sites of non-uniform editing.

Homing guide systems are parameterized similarly to nuclease systems, with elevated probabilities of insertions and deletions at labeled target sites. However, editing of a site within a larger associated editing window does not preclude subsequent editing events within this target window, thus mimicking the maintained compatibility between mutated target sites and gRNA.

In prime editing systems, insertion probabilities are increased with default editing windows. Deletion probabilities are decreased relative to nuclease systems to reflect the increased precision of the single-strand nick mechanism over double stranded breaks. Insertion sequences are drawn from a fixed library.

### Mitochondrial mutation, organelle dynamics, and fitness

Unlike prospective barcodes that begin as a homogeneous population of unedited matrices, mitochondrial genomes exhibit underlying heterogeneity in the founder cell. We allow for a user-provided fraction of mitochondrial genome positions to exhibit some degree of heteroplasmy across all genomes in the founder cell. For simplicity, we assume that mitochondrial variants exist as single-base point mutations. To account for the potential impact of heteroplasmy on cellular fitness, we associate each variant with a severity score drawn from a two-component Gaussian mixture density, with component means μ1 = -1 and μ2 = 1 (arbitrary means corresponding to deleterious and protective scores), a shared standard deviation σ provided by the user that controls the shape of the density curve between μ1 and μ2, and a mixing weight α (also user-provided) representing the fraction of variants that are negatively associated with cellular fitness.

We compute a cell-level heteroplasmy score across k unique position-by-mutation combinations in the genome as the sum, over the k combinations, of the fraction of mitochondrial genomes in the cell carrying mutation i (f_i) multiplied by that mutation's severity score (s_i): S_c = Σ f_i · s_i.

At each cell division point, we allow mitochondrial fusion and fission events to occur. We simplify dynamics such that any two mitochondria within a cell can merge their genomes into one heritable unit, and any one mitochondrion can split its genomes into two mitochondria. Formally, we update our previously defined mitochondria-to-mitochondrial-genome map φ^-1:

During a fusion event, for two fusing mitochondria a and b joining to form mitochondrion c:
- φ^-1(c) = φ^-1(a) union φ^-1(b)
- φ^-1(m) is unchanged for all m other than a, b, and c
- the new mitochondrial count n'_mito = n_mito - 1

During a splitting event of mitochondrion m into genome index sets I_m^1 and I_m^2 (whose union is I_m):
- φ^-1(m1) = I_m^1
- φ^-1(m2) = I_m^2
- φ^-1(r) is unchanged for all r other than m1 and m2
- the new mitochondrial count n'_mito = n_mito + 1
- and φ^-1(m) is removed.

### Mitochondrial inheritance in the native simulator

All mitochondria are replicated before a cell divides. The maintained native simulator then allocates the resulting mitochondrial units stochastically to the daughters and reindexes the daughter-specific mitochondrion-to-genome maps. Configuration values for non-random inheritance are rejected explicitly. Biased, progressive-bottleneck, and clustered inheritance regimes explored during earlier model development are therefore not treated as capabilities of the current implementation; results that depend on those regimes require reimplementation and validation before inclusion in the primary analysis.

### Heteroplasmy-dependent cell fitness

Following this segregation, we transform cell-specific heteroplasmy scores into heteroplasmy-driven cell death probabilities. We formulate the probability that a cell c with heteroplasmy score S_c survives a heteroplasmy checkpoint during cell division using a logistic model: P(survive | S_c) = 1 / (1 + exp[-(α + β(S_c - S_0))]), where S_0 is the heteroplasmy score of the founding cell. The intercept α of the logistic model is obtained from the survival probability of the founder cell, α = log(p0 / (1 - p0)), giving P(survive | S_c) = 1 / (1 + exp[-(log(p0/(1-p0)) + β(S_c - S_0))]); we fix β = 1 such that more positive heteroplasmy scores are associated with increased survival probabilities.

We assume selection by heteroplasmy operates only at the cellular level rather than at the mitochondrial level and do not allow for the selective proliferation of lower-heteroplasmy mitochondrial genomes within a mitochondrion. We do permit the probabilistic loss of some fraction of mitochondria at each cell division point, independent of heteroplasmy.

### Fixed-tree lineage import and recording replay

The fixed-tree workflow accepts the PhysiCell division-event file `cell_lineage.csv`, in which each row records division time, retained parent identifier, and new daughter identifier. Because PhysiCell allows a parent identifier to persist after division, each event is expanded into two new branch segments: a continuation of the retained parent and a branch for the new daughter. Repeating this procedure in chronological order yields an event-resolved binary tree with explicit start time, end time, parent segment, and PhysiCell identifier for every branch. The optional final-cell file `lineage_table.csv` is used to distinguish extant terminal cells from branches removed before the final snapshot. When the final-cell file is omitted, every identifier active at the end of the division log is treated as sampled.

Engineered recording is replayed over the exact duration of every branch. Per-cell-cycle editing probabilities in the simulation configuration are converted to continuous-time hazards, and branch-level event probabilities are obtained from those hazards and elapsed branch time. Edits are irreversible and are inherited by both descendants. If a branch crosses the configured global induction time, it is split into pre-induction and post-induction intervals and the corresponding background or induced hazards are applied to each interval.

Mitochondrial replay uses a separate, deliberately simpler model than the native organelle simulator. Each cell holds a fixed number of sparse mitochondrial genomes. At division, each child resamples that number of genomes with replacement from its parent, representing a fixed-size segregation bottleneck, and inherits new mutations drawn from continuous-time branch hazards. This replay does not model explicit mitochondrial fusion, fission, mitophagy, or variable copy number.

The staged driver copies the lineage-enabled PhysiCell `tumor_3D_lineage` project into a run-specific build directory, compiles it, and terminates growth at a configurable current-cell threshold without modifying the PhysiCell checkout. Intermediate full snapshots are disabled for the large run, while the complete division log and final-cell table are retained. The default production configuration targets 10,000 current cells and then invokes the fixed-tree simulator for both base-editing and mitochondrial modalities.

### Reference-conditioned transcriptome simulation

The optional transcriptome stage uses scDesign3 to fit gene-wise count distributions and gene-gene dependence to a user-provided single-cell reference (Song et al. 2024). References may be supplied as a `SingleCellExperiment` or Seurat RDS object, or as an AnnData H5AD file. The count assay and reference metadata are standardized before fitting. A reference cell-type column is required; lineage pseudotime, two-dimensional spatial coordinates, and additional continuous or categorical covariates may be included when matching fields exist in both the reference and simulated metadata.

For each terminal PhysiCell cell, the workflow constructs a shared metadata record containing the terminal sample identifier, cell type, lineage depth, normalized lineage pseudotime, birth and sampling time, branch length, three-dimensional position, radial distance from the tumor origin, neighbor count, barcode edit count and fraction, and mitochondrial variant and heteroplasmy summaries. Lineage pseudotime is min-max scaled to the range of the reference pseudotime when requested. Similarly, simulated x and y coordinates are scaled to the two supplied reference spatial-coordinate ranges. Other requested covariates are used without automatic transformation and must be represented in both datasets.

The standard covariate-linked implementation calls scDesign3's `construct_data`, `fit_marginal`, and `fit_copula` procedures using negative-binomial marginals and a Gaussian copula, caches the fitted model by reference and formula configuration, and generates new counts with `extract_para` and `simu_new`. The default mean formula includes cell type and optionally a smooth lineage-pseudotime term, a two-dimensional spatial smooth, and user-selected matched covariates. The generated count matrix is stored as a `SingleCellExperiment` whose column names and metadata row names exactly match the terminal identifiers in the lineage-recording outputs. This model links expression to lineage, space, and recorder burden through measured covariates; conditional on those covariates, it does not add a separate phylogenetic covariance process among relatives.

### Molecular recovery, score matrices, and tree reconstruction

At the conclusion of the cell population simulation process, we allow for the selective downsampling of recovered terminal cells to provided cell type-specific recovery probabilities. Then, according to barcode integration and mitochondrial recovery probabilities, we further subset M_bc and M_mt, respectively, for each recovered terminal cell. During the mitochondria recovery process, we use φ_c^-1 to extract the respective genome indices for each recovered mitochondrion. Note that we assume these recovery rates are uniform across integrations and mitochondria regardless of the underlying mutation burden. Thus, we assume missing mutational data is an uninformative feature in downstream tree building steps.

We use a look-up approach to map unique mutations to the cells containing those mutations. For barcode matrices, these unique mutations are identified by collapsing on integration number, genomic position, and mutation identity; for mitochondrial matrices, only genomic position and mutation identity are used to identify unique mutations, and allelic fractions for each variant at each genomic position are calculated. For a population with n_term terminal cells and u_mt and u_bc unique barcode and mitochondrial mutation combinations, we define two score matrices:

- W_mt: an n_term-by-u_mt matrix with entries in [0, 1] (continuous allelic-fraction scores for the mitochondrial mutations), and
- W_bc: an n_term-by-u_bc matrix with binary (0/1) entries for the barcode mutations.

For a given allelic fraction threshold f_af, we further process the mitochondrial score matrix into a binary matrix by setting entry (i,j) to 1 if the corresponding allelic fraction A_ij is at least f_af, and 0 otherwise.

To jointly infer lineage from barcode and mitochondrial modalities, the respective mutation matrices of the paired modalities are concatenated to one another along a cell axis. These modality-specific and joint binary matrices serve as the input to IQ-TREE for phylogenetic inference.

### Clade definitions and overlap metrics

Let a tree be denoted by its vertex set V and edge set E, rooted at v_root. Let V_int be the set of internal (non-leaf) nodes and V_leaf be the set of leaf nodes.

Three approaches are used to map each leaf to a clade. In the topological approach, the tree is cut into a predetermined number of subtrees, and each leaf v is assigned to a clade according to its nearest ancestral cut node. In the mutational signature approach, a fixed number of clades are defined by predominant mutational signatures, and rarer signatures are subsequently merged into the nearest existing clade according to Hamming distance of precomputed score matrices. In the random approach, leaves are mapped randomly to any one of a predefined set of clades.

### Expectation-maximization classification of cell-state transition regimes

To illustrate the information capacity of mitochondria-derived clade assignments, we develop a tree model with observed node states and latent transition process assignments. Unlike a standard hidden Markov tree model, the cell states themselves are treated as observed at the leaves, and the latent variables instead determine which of two transition candidate matrices, Phi0 and Phi1 (uninduced and induced, respectively), generated each internal subtree.

For each internal node v, define subtree X_v as the subtree rooted at v and containing all descendants of v, and introduce a latent binary variable Z_v, where Z_v = 0 if X_v evolved under Phi0, and Z_v = 1 if X_v evolved under Phi1.

Each clade c has a mixture parameter p_c = P(Z_u = 1 | u in c) that reflects the probability that the divisions in clade c use the induced transition matrix Phi1. The EM procedure iteratively seeks to estimate each p_c while keeping Phi0 and Phi1 fixed.

We assume that each leaf node v has an observed cell state x_k (one of K possible states). We probabilistically encode these leaf cell types by defining a one-hot likelihood vector for each leaf v: the i-th entry L_v^i equals 1 if i = x_k, and 0 otherwise.

For an internal node u with child v, the child likelihood vector is propagated upward through each transition matrix. For process r (0 or 1), the propagated likelihood M_r^v(i) = P(X_v | z_u = i, Phi_r), computed as the matrix product of Phi_r (transposed) and the likelihood vector L_v.

We assume independence between children given a parent state. Thus the subtree likelihood for process r at node u is the product, over all children v of u, of M_r^u(i). After log-transforming, this becomes the sum over children of log M_r^u(i).

To prevent numerical underflow, we rescale so that exponentiated log-likelihoods lie on [0, 1], by subtracting the maximum likelihood element of the uninduced transition matrix's log-likelihood vector from each entry.

Since subtrees rooted at each internal node may give rise to child nodes in multiple clades, we define a clade-weighted induced prior p_v for each internal node v according to the relative enrichment of cell types across all downstream clades, representing P(Z_v = 1). The subtree likelihood at node v is then modeled as a weighted mixture of the uninduced and induced likelihoods: P(X_v | z_v = i) = (1 - p_v)·L_0^v(i) + p_v·L_1^v(i) = L_v(i). This likelihood vector is normalized by dividing by its maximum entry to obtain a normalized likelihood.

We fix the root cell state to be x_root,k and define its prior π(i) to equal 1 when i = x_root,k, and 0 otherwise.

The total log-likelihood of all observed leaf states in the tree, log P(X), is given by the log of the sum over i of π(i)·L_root(i), plus the sum over all internal nodes u of a scaling correction term (logfactor_u) determined from the rescaling described above.

For each internal node v, EM finds the posterior probability γ_v = P(Z_v = 1 | X). Given a parent prior vector q, the conditional likelihood under each induction scheme r is P(X | Z_v = r) = the sum over i of q(i)·L_r^v(i), or in log-space, A_r = log of that same sum (computed via the log-sum-exp of log L_r^v(i) weighted by q(i)).

Thus the posterior responsibility that node v contributes to determining whether the induced transition matrix Phi1 was used is γ_v = (p_v·A1) / (p_v·A1 + (1-p_v)·A0).

Each node v is further assigned a tunable weight w_v to determine how strongly it contributes to the M-step. By default, nodes contribute uniformly. 

For each clade c, we accumulate these weighted responsibilities and normalize by the total weight of all nodes that map to c: N_c = sum over v in c of w_v·γ_v, and D_c = sum over v in c of w_v.

The induced fraction for each clade is then updated using a beta prior: p_c = (N_c + α_c - 1) / (D_c + α_c + β_c - 2), which is analogous to the MAP estimate of the beta distribution. The model iterates until convergence and ultimately outputs clade-level induced probabilities.

We threshold clade-specific induction probabilities at 0.5 and label all cells v in clade c with the induction regime y_v: "induced" if p_c ≥ 0.5, and "uninduced" if p_c < 0.5.

We apply this approach to two distinct settings: a cell population that develops with two transition matrix regimes from the beginning, and detection of a state-switching mechanism in a diseased-system context, modeling hematopoietic stem cell development. 

### Sensitivity analysis of input transition matrices

The induction inference workflow requires a priori knowledge of transition matrices for both induced and uninduced conditions, yet this information may not be exactly known in a real world context implementation. To determine the effects of uncertainty in the input transition matrices on EM performance, we add increasing noise to the induced transition matrix according to a mixture parameter 𝛼, such that Phi1_noisy = (1-𝛼)Phi0 + (𝛼)Phi1_true. This sensitivity test captures how the prior passed to the EM modulates model performance as 𝛼 varies from 0 (fully uninformative: every row equal to Phi0) to 1 (fully informative: every row equal to Phi1). 


<!--
Archived pre-restructure Results draft (retained for provenance; not rendered).

We first benchmark TOOL by simulating two prospective lineage tracing systems: a Cas9 nuclease and adenine base editor, parameterized by observed edit counts in FLARE and BASELINE, respectively. Three replicates of each experiment were run, and accuracy scores were mean-aggregated. 



Generally, the base editing system outperformed the nuclease system according to tree reconstruction accuracy, here defined as 1- normalized Robinson-Foulds distance between prospective barcode-estimated trees and ground truth lineage structure. At initial founder integrations above 20 and perfect integration recovery probabilities, BASELINE was able to perfectly reconstruct underlying tree structure (norm RF = 0) across all cell sampling fractions tested once the cell population had reached sufficient size (here, 128 cells). While also improving in accuracy with increasing cell sampling fraction and number of founder integrations, FLARE was not able to reliably capture this perfect tree structure. 

Why does accuracy decrease sometimes with increasing cell population size in baseline


After talking over how we should have a limited number of models to compare, I realized I had a couple of issues with the models from both the 126 run and the 51 run. So, I did a new “12” run, which has four parameter combinations that are run three times apiece. 


We next use TOOL to compare the relative information capacity of mitochondrial lineage tracing systems across a wide parameter space, with the goal of identifying specific parameterizations that perform well on downstream clade reconstruction and classification tasks. 

For simplicity, we seed each simulated cell population with a single founder cell containing 50 mitochondria and an average of five mitochondrial genomes per mitochondrion, each of which has length 16.6 kb. The initial state of mitochondrial heteroplasmy in the founder cell is determined by first parameterizing the fraction of all sites in the mitochondrial genome that can be variable (here, 0.0006). For each site that is designated as a possible heteroplasmy position, the variant allele frequency (VAF) of this variant is drawn from Beta(0.25, 1).


Mitochondria are partitioned into daughter cells according to one of four cell type-specific mitochondrial inheritance mechanisms: 

- Random: daughter 1 inherits a random binomial draw from the parent mitochondrial pool without replacement; daughter 2 inherits the remaining parental mitochondria
- Biased: similar to random, but daughter 1 preferentially receives a larger share of parental mitochondrial pool, as governed by a bias parameter p ≠ 0.5. Daughter cell mitochondrial populations expand or constrict back to the parental number of mitochondria prior to downstream division events.
- Bottleneck: the parent mitochondrial pool is reduced to k founder mitochondria before being split between daughters, such that smaller values of k correspond to larger drift per division. Daughter cell mitochondria pools expand to parental number of mitochondria prior to downstream division events.
- Clustered: mitochondria segregate as fused network clusters of fixed size rather than individually. Following a division event in which each daughter receives half of the parental clusters, the daughter clusters iteratively split according to size until the parental number of network clusters is restored.

Fusion and fission events between mitochondria are encoded by separate rate parameters that govern the expected number of respective events per mitochondrion per cell division. For simplicity, mitochondrial genomes can also be encoded as well-mixed such that, prior to each cell division, all genomes are randomly shuffled across all mitochondria within a cell, thus decoupling variant signals from the unit of inheritance (the mitochondrion). **Well-mixed fusion/fission dynamics are not included in these analyses.

To illustrate the effects on lineage capacity of varying mitochondrial inheritance patterns and dynamics, we used TOOL to simulate four mitochondrial parameter combinations, as shown in Table ______.


All four combinations share the same fusion/fission-rate gradient, which reflects the tendency of more differentiated cell types to experience higher levels of mitochondrial fusion than fission. This gradient is shown in Table ______. 



Model L overrides these parameters for LT and ST (see below). 

Model I is a baseline null model. Each cell type uses a random inheritance pattern such that, after fusion and fission dynamics and doubling, mitochondria are split according to a binomial draw between each daughter cell. Mitochondria counts are permitted to inflate or deflate with time. 

Model J encodes a gradual shift from asymmetric mitochondrial inheritance at early progenitor states to more balanced patterns at differentiated states. Differentiated myeloid and lymphoid cells follow a balanced inheritance pattern that is subject to mitochondrial expansion and contraction to preserve constant mitochondrial counts per cell.

Model K simulates an extreme case of successively weakening bottlenecking events. At each bottleneck, each daughter inherits k mitochondria that are randomly sampled without replacement from the doubled parent pool and subsequently allowed to expand back to the starting mitochondria count. The severe bottleneck at LT guarantees a sharp restriction of genetic diversity across all mitochondria. 

Model L aggregates mitochondria into heritable networks of predefined size at the LT and ST states. During division, half of all networks are passed to each daughter cell such that the original number of mitochondrial clusters is halved. Immediately following this division, the inherited clusters in the daughter cells again split to restore the initial number of clusters for that state. Fusion and fission rates are encoded to be lower in LT and ST to further underscore the heritability of these mitochondrial clusters. Differentiated states (after progenitors) share parameters with Model J. 

Each model was run with three replicates to capture some of the stochasticity that exists in the underlying mutation and mitochondrial inheritance processes. 

Lineage trees reconstructed using only mitochondrial mutational data for each timepoint are plotted in Figure _____. Clades are annotated from innermost to outermost layer according to differentiation induction status, cell type, five-clade GTTC, and MMPC. Because the initial population was seeded with both induced and uninduced progenitor cells at the first cell division, roughly half of all cells in each population arise from the differentiated and undifferentiated transition matrix dynamics. 





The number of callable mitochondrial mutations per cell at varying timepoints and VAF thresholds varies across models, as shown in Figure ____. The unconstrained random inheritance patterns of Model J enables the calling of ~10 mutations per cell with zero thresholding, of which ~4 are highly penetrant and remain after imposing a strict 0.1 VAF threshold. This reduction differs from the progressive bottleneck of Model J, which maintains similarly low counts of highly penetrant variants across all VAF thresholds and timepoints. The unique binary variant patterns across all cells and variants for each model at a late-simulation timepoint are shown in Figure ____. 





Mitochondrial lineage tracing efforts generally seek to identify clonal relationships between cells. We used averaged F1 purity and ARI (together, “overlap metrics”) to quantify the degree of overlap between MMPCs assigned using mitochondrial mutational profiles and GTTCs assigned using ground truth tree topology structure. We manually specify the fixed number of clades in the GTTC approach as 5 and 20 so that MPPC accuracy can be benchmarked against a coarse-grained clonal resolution and a finer-grained subclonal resolution.



The distribution of the maximal number of distinct detectable mitochondrial mutation sites (assuming perfect genome recovery and zero VAF thresholding) across each of the model frameworks is plotted in Figure ___. Model K consistently expresses fewer variant sites because of the severe bottlenecking at early progenitor timepoints. 

The progressive bottlenecking approach of Model K consistently outperformed all other models across both overlap metrics and both GTTC resolutions. The 5-clade and 20-clade GTTC ARI values for Models I, J, and L are strongly correlated with one another, suggesting a consistent albeit weak ability to recapitulate ground truth clade assignments across a range of tree resolutions. Meanwhile, Model K tended to exhibit higher ARI in the 5-clade GTTC setting than in 20-clade GTTC, suggesting that the lower callable mitochondrial mutation count yielded a preferential advantage toward coarse-grained clade assignments over fine-tuned clade placements. This result is also partially explained by the fewer number of MMPCs in Model K than in the three remaining models. The 5-clade GTTC ARI (0.224 +/- 0.085 SEM) is roughly ten times higher than that of the next highest ARI across Models I, J, and L (Model L: 0.025 +/- 0.003 SEM). The overall magnitude of the overlap metrics across Models I, J, and L are similarly underpowered across 5-clade and 20-clade GTTC conditions.

F1 purity scores decrease from 5-clade to 20-clade GTTC across all model iterations, reflecting the increased difficulty of assigning cells to specific subclades using only mitochondrial mutation profiles. 

Smaller numbers of MMPCs in Model K appear to be associated with higher overlap metric scores, whereas no such relationship appears to exist in any of the other three models. 

Cell type-specific mitochondrial inheritance patterns along the differentiation axis enable variable clade reconstruction accuracy across timepoints, as shown in Figure ____. All models generally improve their clade reconstruction potential over time, but the variance in the bottleneck of Model K is much higher across all timepoints; this is concordant with the high standard error observed above, particularly in the 5-clade GTTC setting.

We next sought to determine whether the lineage signal present in mitochondrially-defined lineage trees could perform comparably to ground truth signal in a contrived clade-level classification problem. The expectation maximization workflow defined in ____ uses relative cell type abundances in MMPCs and GTTCs, as well as the underlying mitochondrial tree structure, to estimate, for each defined clade, which of two transition matrices (induced or uninduced, for simplicity) was more likely to have given rise to the set of cell types in that respective clade. Ultimately this model is a cell type enrichment approach grounded in tree structure that relies on the co-occurrence of cell types in the same clades. It raises the testable hypothesis of whether mitochondrial mutation profiles can extend beyond the clonal resolution and to the level of single cell lineage tree reconstruction. 



Because the underlying differences between induced and uninduced matrices lie only in the directionality of the downstream fate bias of intermediate progenitors, relative cell type proportions remain approximately equal between the two regimes across all timepoints (Figure ____). Traditional cell type abundance workflows would thus be underpowered for this scenario.

** I found a bug in the 20-clade GTTC results while I was writing this, so I am only reporting 5-clade GTTC results here for now (despite their being in the figures – I’ll add comparisons between 5-clade and 20-clade once the EM is re-run on 20-clade, then update the figures)

Given that cell type composition within and between clades is subject to random walks along the induced and uninduced transition matrices, we benchmark model performance of MMPCs against that of GTTCs. For each considered EM metric x (namely, ROC-AUC, PR-AUC, F1 score, and accuracy), we compute MMPC_x / GTTC_x.   Across most EM metrics, the progressive bottlenecking approach of Model K yields the highest AUROC relative to the 5-clade GTTC (0.656 +/- 0.04) as well has the highest raw ROC-AUC (0.593 +/- 0.021). PR-AUC is the only relative EM metric on which Model K (0.743 +/- 0.034) does not score best, as the biased inheritance pattern Model J reaches a PR-AUC of 0.815 +/- 0.143. 

** a weak correlation might exist between ari and relative AUROC metric? Need to include more data points to find out

Since the EM assumes a priori knowledge of underlying transition matrix values, we iteratively added a smoothing factor to the induced matrix such that, at alpha = 0, it is identical to the uninduced matrix; at alpha = 0.5, it is an equal-weighted blend of the two; at alpha = 1, it retains its original parameterization. We see that in both GTTC- and MPPC-defined clades, all alpha values greater than zero yield roughly the same ROC-AUC scores within each model. The delta ROC-AUC between alpha = 0 and alpha = 0.1 (the smallest tested alpha value) was greater than the delta between alpha = 0.1 and alpha = 1 in all but Model K. Taken together, these results suggest that the EM model is robust to uncertainties in exact transition matrix parameterizations between the two regimes, so long as the discriminative signal between them points toward the correct underlying biology. 

The low magnitude of each EM metric, together with the low relative performance against GTTC gold standard, suggest that the tree structures inferred by mitochondrial mutation profiles alone are not sufficiently powered for downstream inference.


End results section writeup 7/22

# Results

The tree reconstruction accuracy A_ℓ for tree ℓ is defined as A_ℓ = 1 - RF_norm, where RF_norm is the normalized Robinson-Foulds distance between ℓ and the known ground truth phylogeny, subsetted to only include cells whose profiles were recovered in the post-processing steps and included in ℓ.

Add mito fallelic fraction thresholding of only random inheritance patterns here. 


**** fix the colors in the alpha sweep plot to match the regime color code from earlier


We first measure the reconstructive power of a base editing system across parameters determined at the time of experimental design as well as during experimental workflows. This system assumes each target is defined with high specificity without an adjacent editing window. Figure [XXXX] shows the interdependent effects that the number of barcode integrations in the founder cell, the probability of recovering a given integration in sequencing runs, and the number of base editing targets have on tree reconstructive power over simulation time. Final cell population sizes scaled to approximately 2^12 = 2048 cells. At high integration counts, the base editing system achieves near- or perfect accuracy across even non-1 integration recovery fractions. Critically, the number of initial founder integrations accounts for much of the variance in accuracy scores.

TODO:
We next explored whether sufficient signal exists in mitochondrial mutational profiles to improve tree reconstruction accuracy in cases where experimental circumstances limit the power of prospective barcoding approaches alone. We assume a mitochondrial capture rate of 50% and generate binary score matrix inputs to IQ-Tree by thresholding mitochondrial signatures at allelic fractions of 0, 0.1, and 0.2. Independently, these mutational profiles fail to capture much of the ground truth tree topology.



Induced vs uninduced clade inference

We first define seven unique mitochondrial inheritance regimes, detailed in Table _____. These regimes were chosen to capture a diversity of selective pressures across varying levels of model complexity.



In all, 42 unique simulation parameterizations were generated by iterating through each unique combination of mitochondrial heteroplasmy, fusion and fission rate, and inheritance model regimes. Each parameter combination is run three times. Results of all 42 runs are shown in Supplementary Figure ____. 

We first assessed whether clades defined by binary single cell mitochondrial allele matrices recapitulated those defined by known tree topology. For each simulation run, every terminal leaf was assigned to a mitochondrially-defined mutation profile clade (MMPC) by grouping together cells with identical mitochondrial haplotypes. If a haplotype is shared by at least five cells, it defines its own clades; otherwise, cells with a rarer haplotype are absorbed into whichever sufficiently large clade their mutational profiles resemble. MMPCs are agnostic to reconstructed tree structure.

Ground truth topological clades (GTTCs) defined using the topology of the underlying ground truth tree were generated by iteratively splitting the largest remaining ground truth sublineage at an internal cut-node until a predefined number of N subtrees remain; here, GTTCs were defined for N = 5 and N = 20 to capture varying degrees of phylogenetic resolution. Each terminal leaf was assigned to the clade of the earliest ancestor that falls at or below a cut point. With two distinct N values, we can better determine whether certain parameter combinations more effectively capture global topological variation or more subtle subclonal differences.


An aggregated purity metric was used to quantify the strength of association between a given simulation run’s MMPCs and GTTCs. For each GTTC clade k, define n(k, j) as the number of cells in k that fall into each MMPC clade j. The dominant MMPC clade j for GTTC clade k is defined as max_j[n(k, j)]. GTTC → MMPC purity is defined by aggregating across all k in K such that:

Purity_{GTTC → MMPC} = sum(over k in K) [[ max_j[n(k,j)] / N]] 

where N is the total number of cells.

Similarly, we find MMPC → GTTC purity:

Purity_{MMPC → GTTC} = sum(over j in J) [[ max_k[n(j,k)] / N]] 

To penalize mappings that result in low purity in one direction and high in the other, we define aggregated purity between a run’s MMPC and GTTC sets as the harmonic mean F of Purity_{GTTC → MMPC}, Purity_{MMPC → GTTC}: 

F1 = 2*Purity_{GTTC → MMPC}*Purity_{MMPC → GTTC} / (Purity_{Purity_{GTTC → MMPC} + MMPC → GTTC})

Comparing Purity_{MMPC → GTTC} and Purity_{GTTC → MMPC} values for a given parameter combination has the added benefit of revealing whether terminal leaves are over- or under-clustered. When Purity_{MMPC} > Purity_{GTTC}, the MMPC set has produced too many clades, and each MMPC is mostly composed of cells from a single GTTC, but each GTTC is split across multiple MMPCs (and vice versa). Generally, this under-clustering is the dominant pattern observed across the tested inheritance regimes, but **CERTAIN REGIMES HAVE A HIGHER RATE OF UNDERCLUSTERING – EXPLAIN. See Figure 3. 

We marginalized F1 purity scores across mitochondrial inheritance regimes, rates of fusion and fission events, and whether mitochondrial heteroplasmy variants are associated with a small fitness penalty at the coarse-level GT-5 clade level and the more fine-grained GT-20 clade level.

We sought to benchmark the extent to which MMPCs could recapitulate results from GTTCs. Due to the stochasticity of the cell type heterogeneity within each topological clade (according to predefined transition matrices), we defined, for each parameter combination, a score of clade-wise induction prediction accuracy that is relative to ground truth prediction accuracy as X_MMPC/X_GTTC, where X is a model evaluation metric (typically AU-ROC). By definition, GTTCs perfectly reflect underlying tree topology, and they were shown to perform well on the EM classification task. 

We deploy the EM inference pipeline to a model of HSC differentiation that has six discrete cell states governed by transition matrix dynamics shown in Figure _____. Briefly, the simulator begins with a founder cell in the long-term hematopoietic stem cell state (LT), followed by the short-term hematopoietic stem cell (ST), and a branch point toward myeloid multipotent progenitor (MyMPP) and terminal myeloid (Mye) cells or toward multipotent lymphoid progenitors (LyMPP) and terminal lymphoid cells. Exact transition matrix values are provided in Figure 4. 

Coefficient of variation was used to quantify the stability of F1 purity and AUROC estimates across parameter combinations (DID NOT FIT IN FIGURE, MAYBE REMOVE). 

We elected to distill the 42 parameter combination space into a five parameter combination subset that 

We next looked to identify the relationship between mitochondrial clade organization (as defined by the F1) and performance on the EM benchmarking task. WEAK ASSOCIATION, FIND CORRELATION BETWEEN ALL 42 COMBOS (TODO)




The sensitivity of the EM classification scheme to input transition matrix values was computed by varying the induced matrix as a linear combination of the induced and uninduced values. Since Phi0 and Phi1 differ only in their pairwise [Myeloid Progenitor, Lymphoid Progenitor] → [Myeloid, Lymphoid] transition probabilities, the scaling only affects these respective values. Resulting Phi0 and Phi1 matrices are shown in Figure ______. When input induced and uninduced transition matrices are identical at 𝛼 = 0, the model effectively picks at random; when a minimal fraction of the true directional difference is injected at 𝛼 = 0.25, performance nearly reaches its ceiling. This suggests that the EM is sensitive to weak directional contrasts between input induced and uninduced matrices and does not require high-accuracy transition estimates to accurately discriminate between transition matrix regimes. 


-->

## Discussion

### Separating lineage generation, molecular recording, and observation

This framework recapitulates the effects of experimental and computational limitations on single-cell lineage reconstruction while placing prospective and retrospective recording modalities in the same simulated system. Its main architectural contribution is the separation of the lineage-generating process from the molecular record written on that lineage and the terminal observation generated from each sampled cell. The native simulator remains useful for controlled experiments in which population state, mitochondrial dynamics, and recorder activity must be varied jointly. Fixed-tree replay extends the same molecular experiments to lineages generated by another model, demonstrated here with spatial tumor growth in PhysiCell. The scDesign3 stage then provides a reference-conditioned expression phenotype aligned to the same terminal-cell identifiers.

This separation also clarifies which comparison is being made. Differences between native simulations and PhysiCell replay can arise from lineage topology and branch-time distributions, whereas recorder parameter sweeps on a single imported tree isolate the molecular observation process. Reusing a fixed PhysiCell tree across recorder configurations should therefore provide a controlled way to evaluate how edit rate, induction timing, molecular recovery, and mitochondrial bottleneck size affect reconstruction in a spatially generated tumor lineage.

### Model assumptions

Mutational processes, stochastic and heteroplasmy-driven cell-death events, and cell-type transitions in the native simulator are parameterized in a discrete framework in which events occur at the temporal resolution Δt. With sufficiently small Δt, these processes approximate continuous-time behavior. If intermediate data are requested strictly between nΔt and (n+1)Δt, mutation profile, death status, and cell type reflect the state available at nΔt. In contrast, fixed-tree replay converts per-cycle probabilities to continuous-time hazards and applies them over each imported branch duration. These distinct time formulations should not be interpreted as identical generative models.

The explicit parameterization of mutation, population, and heteroplasmy dynamics provides direct control over simulated lineage-recording experiments. Where possible, these probability distributions should be calibrated against empirical data from the relevant cell line, recorder construct, and sequencing protocol. The current native mitochondrial implementation permits stochastic fusion and fission but neglects higher-order network structure and uses random daughter allocation. The fixed-tree mitochondrial replay is more abstract still: it represents inheritance as resampling a constant number of genomes at every division and does not model organelles explicitly.

Changes in cellular fitness conferred by mitochondrial variants are estimated from variant allele fractions across all mitochondrial genomes in a cell rather than from gross mutational load. This normalization buffers sudden changes in fitness score that might otherwise arise from random allele segregation. Selection operates only at the cell level, and the model does not represent within-cell competition among mitochondrial genomes.

### PhysiCell and transcriptomic extensions

PhysiCell supplies three-dimensional tumor growth, mechanical interactions, and irregular division timing that are absent from the native lineage generator (Ghaffarizadeh et al. 2018). The current adapter consumes the complete division history and the final current-cell table, but only a subset of terminal spatial metadata is exported into downstream analyses. Adding phenotype, microenvironmental substrate concentrations, cell-cycle phase, and death events would allow recording and expression to depend directly on tumor ecology rather than on lineage depth and position alone.

The scDesign3 integration links simulated expression to covariates shared with a real reference, including cell type, lineage pseudotime, space, and recording-derived burden (Song et al. 2024). It does not infer a biological effect for a covariate absent from the reference, and a recording statistic influences expression only when that statistic or a defensible proxy is available in the reference metadata and included in the model formula. The standard implementation also samples terminal expression conditionally on those covariates without an additional tree-indexed residual process. Consequently, related cells can resemble one another through inherited or spatial covariates, but not solely because they share recent ancestry. A future phylogenetic residual model could introduce that dependence explicitly.

### Limitations and next analyses

Several Results subsections remain provisional. The 10,000-cell PhysiCell workflow is configured but has not yet been used for the biological benchmark reported here, and the current validation establishes interoperability rather than large-scale biological realism. Historical native analyses involving biased, progressive-bottleneck, or clustered mitochondrial inheritance require a restored and tested implementation before they can support primary claims. The known 20-clade analysis error must also be rerun. Finally, reference-conditioned expression should be evaluated with held-out reference cells and lineage-aware summary statistics before being used to benchmark downstream transcriptomic methods.

Together, these limitations define a practical next analysis: generate replicate 10,000-cell PhysiCell tumors, replay paired engineered and mitochondrial recorders across the same trees, synthesize covariate-linked transcriptomes from a biologically matched tumor reference, and quantify reconstruction and cell-state recovery while varying only one simulation layer at a time.

### Code availability

The implementation is written primarily in R, with shell drivers for staged PhysiCell execution, and is available at [INSERT GITHUB LINK]. Configuration files and command-line workflows required to reproduce the native, PhysiCell replay, and optional scDesign3 stages are included with the source.

## References

Anzalone, Andrew V., Luke W. Koblan, and David R. Liu. 2020. "Genome Editing with CRISPR-Cas Nucleases, Base Editors, Transposases and Prime Editors." Nature Biotechnology 38 (7): 824-44.

Anzalone, Andrew V., Peyton B. Randolph, Jessie R. Davis, Alexander A. Sousa, Luke W. Koblan, Jonathan M. Levy, Peter J. Chen, et al. 2019. "Search-and-Replace Genome Editing without Double-Strand Breaks or Donor DNA." Nature 576 (7785): 149-57.

Chen, Peter J., and David R. Liu. 2023. "Prime Editing for Precise and Highly Versatile Genome Manipulation." Nature Reviews. Genetics 24 (3): 161-77.

Chen, Wen, Huakan Zhao, and Yongsheng Li. 2023. "Mitochondrial Dynamics in Health and Disease: Mechanisms and Potential Targets." Signal Transduction and Targeted Therapy 8 (1): 333.

Elson, J. L., D. C. Samuels, D. M. Turnbull, and P. F. Chinnery. 2001. "Random Intracellular Drift Explains the Clonal Expansion of Mitochondrial DNA Mutations with Age." American Journal of Human Genetics 68 (3): 802-6.

Frieda, Kirsten L., James M. Linton, Sahand Hormoz, Joonhyuk Choi, Ke-Huan K. Chow, Zakary S. Singer, Mark W. Budde, Michael B. Elowitz, and Long Cai. 2017. "Synthetic Recording and in Situ Readout of Lineage Information in Single Cells." Nature 541 (7635): 107-11.

Ghaffarizadeh, Ahmadreza, Randy Heiland, Samuel H. Friedman, Shannon M. Mumenthaler, and Paul Macklin. 2018. "PhysiCell: An Open Source Physics-Based Cell Simulator for 3-D Multicellular Systems." PLOS Computational Biology 14 (2): e1005991. https://doi.org/10.1371/journal.pcbi.1005991.

Hirose, Misa, Paul Schilf, Yask Gupta, Kim Zarse, Axel Kunstner, Anke Fahnrich, Hauke Busch, et al. 2018. "Low-Level Mitochondrial Heteroplasmy Modulates DNA Replication, Glucose Metabolism and Lifespan in Mice." Scientific Reports 8 (1): 5872.

Jones, Matthew G., Alex Khodaverdian, Jeffrey J. Quinn, Michelle M. Chan, Jeffrey A. Hussmann, Robert Wang, Chenling Xu, Jonathan S. Weissman, and Nir Yosef. 2020. "Inference of Single-Cell Phylogenies from Lineage Tracing Data Using Cassiopeia." Genome Biology 21 (1): 92.

Kalhor, Reza, Kian Kalhor, Leo Mejia, Kathleen Leeper, Amanda Graveline, Prashant Mali, and George M. Church. 2018. "Developmental Barcoding of Whole Mouse via Homing CRISPR." Science (New York, N.Y.) 361 (6405): eaat9804.

Kang, Eunju, Xinjian Wang, Rebecca Tippner-Hedges, Hong Ma, Clifford D. L. Folmes, Nuria Marti Gutierrez, Yeonmi Lee, et al. 2016. "Age-Related Accumulation of Somatic Mitochondrial DNA Mutations in Adult-Derived Human IPSCs." Cell Stem Cell 18 (5): 625-36.

Komor, Alexis C., Yongjoo B. Kim, Michael S. Packer, John A. Zuris, and David R. Liu. 2016. "Programmable Editing of a Target Base in Genomic DNA without Double-Stranded DNA Cleavage." Nature 533 (7603): 420-24.

Kwok, Aaron Wing Cheung, Chen Qiao, Rongting Huang, Mai-Har Sham, Joshua W. K. Ho, and Yuanhua Huang. 2022. "MQuad Enables Clonal Substructure Discovery Using Single Cell Mitochondrial Variants." Nature Communications 13 (1): 1205.

Lareau, Caleb A., Leif S. Ludwig, Christoph Muus, Satyen H. Gohil, Tongtong Zhao, Zachary Chiang, Karin Pelka, et al. 2021. "Massively Parallel Single-Cell Mitochondrial DNA Genotyping and Chromatin Profiling." Nature Biotechnology 39 (4): 451-61.

Li, Siqi, Kun Wang, Xin Wang, and Zheng Hu. 2026. "Single-Cell Mitochondrial Lineage Tracing: Opportunities and Challenges." Quantitative Biology (Beijing, China) 14 (1): e70018.

Liu, Fengshuo, Xiang Zhang, and Yipeng Yang. 2024. "Simulation of CRISPR-Cas9 Editing on Evolving Barcode and Accuracy of Lineage Tracing." Scientific Reports 14 (1): 19213.

Ludwig, Leif S., Caleb A. Lareau, Jacob C. Ulirsch, Elena Christian, Christoph Muus, Lauren H. Li, Karin Pelka, et al. 2019. "Lineage Tracing in Humans Enabled by Mitochondrial Mutations and Single-Cell Genomics." Cell. https://doi.org/10.1016/j.cell.2019.01.022.

Mao, Shanjun, Chenyang Zhang, Runjiu Chen, Shan Tang, Xiaodan Fan, and Jie Hu. 2025. "Cell Lineage Tracing: Methods, Applications, and Challenges." Quantitative Biology (Beijing, China) 13 (4): e70006.

McKenna, Aaron, Gregory M. Findlay, James A. Gagnon, Marshall S. Horwitz, Alexander F. Schier, and Jay Shendure. 2016. "Whole-Organism Lineage Tracing by Combinatorial and Cumulative Genome Editing." Science 353 (6298): aaf7907.

McKenna, Aaron, and James A. Gagnon. 2019. "Recording Development with Single Cell Dynamic Lineage Tracing." Development (Cambridge, England) 146 (12): dev169730.

Porto, Elizabeth M., Alexis C. Komor, Ian M. Slaymaker, and Gene W. Yeo. 2020. "Base Editing: Advances and Therapeutic Opportunities." Nature Reviews. Drug Discovery 19 (12): 839-59.

Salvador-Martinez, Irepan, Marco Grillo, Michalis Averof, and Maximilian J. Telford. 2019. "Is It Possible to Reconstruct an Accurate Cell Lineage Using CRISPR Recorders?" ELife 8 (January). https://doi.org/10.7554/eLife.40292.

Song, Dongyuan, Qingyang Wang, Guanao Yan, Tianyang Liu, Tianyi Sun, and Jingyi Jessica Li. 2024. "scDesign3 Generates Realistic In Silico Data for Multimodal Single-Cell and Spatial Omics." Nature Biotechnology 42: 247-52. https://doi.org/10.1038/s41587-023-01772-1.

Spanjaard, Bastiaan, Bo Hu, Nina Mitic, Pedro Olivares-Chauvet, Sharan Janjuha, Nikolay Ninov, and Jan Philipp Junker. 2018. "Simultaneous Lineage Tracing and Cell-Type Identification Using CRISPR-Cas9-Induced Genetic Scars." Nature Biotechnology 36 (5): 469-73.

Stewart, James B., and Patrick F. Chinnery. 2015. "The Dynamics of Mitochondrial DNA Heteroplasmy: Implications for Human Health and Disease." Nature Reviews. Genetics 16 (9): 530-42.

Tilokani, Lisa, Shun Nagashima, Vincent Paupe, and Julien Prudent. 2018. "Mitochondrial Dynamics: Overview of Molecular Mechanisms." Essays in Biochemistry 62 (3): 341-60.

Wagner, Daniel E., and Allon M. Klein. 2020. "Lineage Tracing Meets Single-Cell Omics: Opportunities and Challenges." Nature Reviews. Genetics 21 (7): 410-27.

Wallace, Douglas C., and Dimitra Chalkia. 2013. "Mitochondrial DNA Genetics and the Heteroplasmy Conundrum in Evolution and Disease." Cold Spring Harbor Perspectives in Biology 5 (11): a021220.

Weng, Chen, Fulong Yu, Dian Yang, Michael Poeschla, L. Alexander Liggett, Matthew G. Jones, Xiaojie Qiu, et al. 2024. "Deciphering Cell States and Genealogies of Human Haematopoiesis." Nature 627 (8003): 389-98.

Xu, Jin, Kevin Nuno, Ulrike M. Litzenburger, Yanyan Qi, M. Ryan Corces, Ravindra Majeti, and Howard Y. Chang. 2019. "Single-Cell Lineage Tracing by Endogenous Mutations Enriched in Transposase Accessible Mitochondrial DNA." ELife 8 (April). https://doi.org/10.7554/eLife.45105.

Ye, Kaixiong, Jian Lu, Fei Ma, Alon Keinan, and Zhenglong Gu. 2014. "Extensive Pathogenicity of Mitochondrial Heteroplasmy in Healthy Human Individuals." Proceedings of the National Academy of Sciences of the United States of America 111 (29): 10654-59.
