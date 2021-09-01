# CoVEME-tutorial
tutorial adapted from http://beast.community/thorney_beast

## Introduction

The SARS-CoV-2 pandemic has highlighted the importance and challenge of large-scale molecular epidemiology. Through 
Herculean effort, academic and public health labs over the world have undertaken genomic sequencing on a scale that 
dwarfs any previous outbreak response. The magnitude of the resulting dataset makes many standard phylogenetic
approaches impractical. Here, we highlight some recent modifications we have made to make the analysis of tens of thousand of 
taxa in BEAST feasible on a time-scale that can help inform public health responses.

Many of these improvements are based on re-implementations of approaches used in other phylogenetic software, and we
are grateful to the entire phylogenetic community. Incorporating these approaches allows us to leverage the full
suite of models already available in BEAST. While the model used here dates back to Zuckerkandl E, Pauling L. 1962, this 
implementation was inspired by [Didelot et al., 2018](https://doi.org/10.1093/nar/gky783).

<div class="alert alert-success" role="alert"><i class="fa fa-download fa-lg"></i>
All the analyses described here can be performed using the <a href="https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre_thorney_v0.1.2"> BEASTv.1.10.5pre_thorney_v0.1.2 </a>pre-release.
</div>


*NB: To the best of our knowledge Jeff Thorne and colleges work in [multidivtime](https://brcwebportal.cos.ncsu.edu/thorne/multidivtime.html) 
represents the first implementation of related approaches in a Bayesian framework.*

### A note about efficiency 

The cost of an operation in a standard BEAST analysis is the product of how long that operation
takes, and how many times it is called. In any analysis calculating the likelihood of the sequence data over a tree (tree-data-likelihood) 
is the most computationally taxing step. BEAST mitigates this by delegating to the BEAGLE library, which limits this cost by
caching calculations on subtrees that remain unchanged between states, and computing the likelihood in parallel on high 
performance GPUs when available. These approaches reduce the number of times the full likelihood is calculated and
increase the speed of that calculation when it is needed. The analysis can be sped up further by fixing the substitution
model parameters and clock rate, if reasonable estimates exists, as sampling these parameters requires a full 
recalculation.

Although the cost of the tree-data-likelihood can be greatly reduced by fixing the parameters that govern the substitution model, 
estimates of the these parameters are not always transferable between datasets, even datasets from the same
pandemic. This is certainly true of SARS-CoV-2 where estimates of the rate can vary based on the time-scale and context of the 
data. Even when it is valid to fix these parameters, large datasets become challenging for a second reason. 
All GPUs have finite memory, and it is challenging if not impossible to find hardware configurations that 
can cope with alignments of tens of thousands of sequences. 


# Tutorial

Before we begin - this analysis requires a rooted precalculated tree with branchlengths in substitution space. For this tutorial we'll be using a 1,000 tip tree built with iqtree2 from gisaid data up until July 2020.


To start we will use beauti to generate a standard BEAST analysis. We will then adapt this xml to run our "Thorney" analysis.
Although we won't be using an alignment in our likelihood calculations, beauti requires one to load the taxa.
There is a mock alignment in the data directory, which contains a single base for each sequence. 
We will fire up Beauti and import it in the
partitions panel. 

Next we will assign dates to the taxa. In the tips panel select `use tip dates` and then `parse dates`. The dates follow the last `|` in the taxa label.

Because we won't be using we can leave the default HKY site model in the site panel. 

For this analysis we'll be using a strict clock model (although any could be used). Select that in the clocks panel.

We will be using a Skygrid with weekly gridpoints. For this dataset, a skygrid cutoff of 0.75 with 39 gridpoints will do the trick.

Next we'll import our starting tree and selected it as the user specified tree. This can be done using File>Import Data.

For now we'll leave all the priors and operators and export our xml. Next we'll adapt this xml 
to use the new likelihood and data structures.

## An alternative likelihood function

Starting from our standard BEAUTi-generated xml you can delete the `<alignment>` and `<patterns/>` block and instead specify a 
"data tree" and a starting tree that are topologically consistent (every clade in the data tree must be present in the starting tree).

For this analyis we already provided a starting tree with mutations in substitution space. We can copy this tree and call it dataTree as below.
NB: Note we set the rescaled root height equal to 1 year to provide a non Infinite starting likelihood
```xml
	<rescaledTree id="startingTree" height="1.0">
        <newick  usingDates="true">
	        ((taxa1|2020-02-17,...);
        </newick>
    </rescaledTree>

<newick id="dataTree" usingDates="false" usingHeights="true">
	((taxa1|2020-02-17,...);
</newick>
```

The new likelihood function requires a map between the number of mutations along a branch in the data tree and 
the duration of that same branch in the estimated time tree. To ensure this mapping is never broken, we have implemented
a  `<constrainedTreeModel>` that takes the place of the standard `<treeModel>` block. The constrained tree model allows
for sampling different resolutions of polytomies that may be present in the data tree without violating the topological requirements
of the tree likelihood.
```xml
<!-- <treeModel id="treeModel">-->
<!--    <coalescentTree idref="startingTree"/>-->
<!--    <rootHeight>-->
<!--        <parameter id="treeModel.rootHeight"/>-->
<!--    </rootHeight>-->
<!--    <nodeHeights internalNodes="true">-->
<!--        <parameter id="treeModel.internalNodeHeights"/>-->
<!--    </nodeHeights>-->
<!--    <nodeHeights internalNodes="true" rootNode="true">-->
<!--        <parameter id="treeModel.allInternalNodeHeights"/>-->
<!--    </nodeHeights>-->
<!-- </treeModel>-->
<constrainedTreeModel id = "treeModel">
	<tree idref="startingTree"/>
	<constraintsTree>
		<tree idref="dataTree"/>
	</constraintsTree>
</constrainedTreeModel>
        
<!-- Statistic for root height of the tree       -->
<treeHeightStatistic id="treeModel.rootHeight">
    <treeModel idref="treeModel"/>
</treeHeightStatistic>
    --snip--
<nodeHeightOperator type="scaleRoot" weight="3" scaleFactor="0.75" >
    <treeModel idref="treeModel"/>
</nodeHeightOperator>

```
Note that because of how the a `<constrainedTreeModel>` stores its root height we need a
different root height statistic and root height operator than what is used in a standard analysis.



We can now replace the standard `<treeDataLikelihood>` and its components with our new likelihood as indicated by the comments
below.
```xml
<!-- Likelihood for tree given sequence data                                 -->

<!-- The strict clock (Uniform rates across branches)                        -->
<strictClockBranchRates id="branchRates">
    <rate>
        <parameter id="clock.rate" value="0.001"/>
    </rate>
</strictClockBranchRates>

<!--<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">-->
<!--<treeModel idref="treeModel"/>-->
<!--<strictClockBranchRates idref="branchRates"/>-->
<!--</rateStatistic>-->


<!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985)             -->
<!--<HKYModel id="hky">-->
<!--<frequencies>-->
<!-- <frequencyModel dataType="nucleotide">-->
<!--  <frequencies>-->
<!--   <parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>-->
<!--  </frequencies>-->
<!-- </frequencyModel>-->
<!--</frequencies>-->
<!--<kappa>-->
<!-- <parameter id="kappa" value="2.0" lower="0.0"/>-->
<!--</kappa>-->
<!--</HKYModel>-->

<!-- site model                                                              -->
<!--<siteModel id="siteModel">-->
<!--<substitutionModel>-->
<!-- <HKYModel idref="hky"/>-->
<!--</substitutionModel>-->
<!--<gammaShape gammaCategories="4">-->
<!-- <parameter id="alpha" value="0.5" lower="0.0"/>-->
<!--</gammaShape>-->
<!--</siteModel>-->

<!--<statistic id="mu" name="mu">-->
<!--<siteModel idref="siteModel"/>-->
<!--</statistic>-->


<!-- Likelihood for tree given sequence data                                 -->
<!--<treeDataLikelihood id="treeLikelihood" useAmbiguities="false">-->
<!--<partition>-->
<!-- <patterns idref="patterns"/>-->
<!-- <siteModel idref="siteModel"/>-->
<!--</partition>-->
<!--<treeModel idref="treeModel"/>-->
<!--<strictClockBranchRates idref="branchRates"/>-->
<!--</treeDataLikelihood>-->

<thorneyTreeLikelihood id="treeLikelihood">
	<constrainedTreeModel idref="treeModel"/>
	<strictClockBranchRates idref="branchRates"/>
    <branchLengthLikelihood id="branchLengthLikelihood" scale="29903.0"/>
	
	<constrainedBranchLengthProvider scale="29903.0"> 
		<constrainedTreeModel idref="treeModel"/>
		<dataTree>
			<tree idref="dataTree"/>
		</dataTree>
	</constrainedBranchLengthProvider>
</thorneyTreeLikelihood>
```
Here `<branchLengthLikelihood>` accepts any branch rate model -  `<strictClockBranchRates>` in our case. It also scales the substitution
rate from substitutions/site/year to substitutions/year (i.e. `scale` represents the number of sites). 
`<constrainedBranchLengthProvider>` maps the branches in the data tree to those in the time tree and provides the number 
of mutations that occur along that branch. Similar to above, `scale` represents a multiplier to scale the branchlengths
in the data tree. In this example the original alignment used to make the data tree had 29,903 sites.

### A note about topologies

Our reliance on a fixed data tree means that every clade in the data must be present in the time tree. However, if the 
data tree is unresolved and contains polytomies, we can sample from possible resolutions without breaking that mapping.
(In these cases we know the inserted branches have 0 mutations.) The tree operators below accept a `<constrainedTreeModel>` 
and are compatible with these constraints. They allow us to sample node heights and polytomy resolutions within the 
constraints imposed by the data tree.

We can delete the other tree operators in the xml - if present.

*NB: The topological operators below will exit with an error if the constraints tree is fully resolved (i.e. there are
no polytomies over which to sample). In that case only the node height operators can be used.*

```xml
<nodeHeightOperator type="scaleRoot" weight="3" scaleFactor="0.75" >
    <treeModel idref="treeModel"/>
</nodeHeightOperator>           

<nodeHeightOperator type="uniform" weight="30">
    <treeModel idref="treeModel"/>
</nodeHeightOperator>
		
<uniformSubtreePruneRegraft weight="30">
    <treeModel idref="treeModel"/>
</uniformSubtreePruneRegraft>
		
<narrowExchange weight="30">
    <treeModel idref="treeModel"/>
</narrowExchange>
		
<wideExchange weight="10">
    <constratreeModelinedTreeModel idref="treeModel"/>
</wideExchange>
		
<wilsonBalding weight="10">
    <treeModel idref="treeModel"/>
</wilsonBalding>
```
		
## More efficient data structures  

Without the computational cost of the tree-data-likelihood, the routine processes of collating many coalescent
intervals and managing a large tree's internal structure become computational bottlenecks. Although these operations
are relatively efficient they occur at every state in the MCMC chain, and so even small inefficiencies add up.

To increase efficiency on large datasets we have implemented new data structures for the time tree and coalescent 
intervals. These implementations form the basis of the `<constrainedTreeModel>` above, but they can be used in any 
analysis and improve performance in even moderately sized datasets (a few hundred taxa).

The tree class can be used by commenting out the standard treemodel and replacing it as below.

```xml

<!-- Generate a tree model                                                   -->

<!-- <treeModel id="treeModel">-->
<!--    <coalescentTree idref="startingTree"/>-->
<!--    <rootHeight>-->
<!--        <parameter id="treeModel.rootHeight"/>-->
<!--    </rootHeight>-->
<!--    <nodeHeights internalNodes="true">-->
<!--        <parameter id="treeModel.internalNodeHeights"/>-->
<!--    </nodeHeights>-->
<!--    <nodeHeights internalNodes="true" rootNode="true">-->
<!--        <parameter id="treeModel.allInternalNodeHeights"/>-->
<!--    </nodeHeights>-->
<!-- </treeModel>-->

<!-- Generate a tree model                                                   -->
<bigFastTreeModel id="treeModel">
    <tree idref="startingTree"/>
</bigFastTreeModel>
<!-- Statistic for root height of the tree       -->
<treeHeightStatistic id="treeModel.rootHeight">
	<treeModel idref="treeModel"/>
</treeHeightStatistic>
	

--snip--

<!--<scaleOperator scaleFactor="0.75" weight="3">-->
<!--     <parameter idref="treeModel.rootHeight"/>-->
<!-- </scaleOperator>-->
<nodeHeightOperator type="scaleRoot" weight="3" scaleFactor="0.75" >
	<treeModel idref="treeModel"/>
</nodeHeightOperator>
```
Note: As is in the case of the `<constrainedTreeModel>`, this change affects how the root height is stored and so we need a new root height statistic and operator.

 

For our purposes here we will use the new efficient coalescent intervals.
We can add them to any demographic model by replacing the `<populationTree>` in a coalescent model
with `<intervals>`. The example below is for our skygrid coalescent model.

```xml
	<gmrfSkyGridLikelihood id="skygrid">
		<populationSizes>

			<!-- skygrid.logPopSize is in log units unlike other popSize                 -->
			<parameter id="skygrid.logPopSize" dimension="39" value="1.0"/>
		</populationSizes>
		<precisionParameter>
			<parameter id="skygrid.precision" value="0.1" lower="0.0"/>
		</precisionParameter>
		<numGridPoints>
			<parameter id="skygrid.numGridPoints" value="38.0"/>
		</numGridPoints>
		<cutOff>
			<parameter id="skygrid.cutOff" value="0.75"/>
		</cutOff>
		<!-- <populationTree>
			<treeModel idref="treeModel"/>
		</populationTree> -->
        <intervals>
		<bigFastTreeIntervals>
			<treeModel idref="treeModel"/>
		</bigFastTreeIntervals>
	</intervals>
	</gmrfSkyGridLikelihood>
```

Because this is a prerelease and represents in progress work in progress at the momenet we also need
to use a different operator on the skygridPopulation sizes.
```xml

		<!-- <gmrfGridBlockUpdateOperator scaleFactor="1.0" weight="2">
			<gmrfSkyrideLikelihood idref="skygrid"/>
		</gmrfGridBlockUpdateOperator> -->
		<gmrfSkygridBlockUpdateOperator scaleFactor="1.1" weight="10">
			<gmrfSkygridLikelihood idref="skygrid"/>
		</gmrfSkygridBlockUpdateOperator>
```

## Managing large xmls with beastgen

These changes can be made to any xml generated by BEAUti; however, editing large xmls is clunky and tedious. An alternative is to use [BEASTGen](http://beast.community/beastgen), a small command line
utility that generates xmls from an Apache FreeMarker template and an input data file, which can be an alignment, a nexus tree, 
or a BEAST xml. 


To adapt out xml into a flexible template that can be used with a new dataset we need to replace the taxa and tree blocks. The `<taxa>` block can be replaced with 
```xml
	<!-- The list of taxa to be analysed (can also include dates/ages).          -->
	<!-- ntax=${taxa?size}                                               -->
	<taxa id="taxa">
<#list taxa as taxon>
		<taxon id="${taxon.id}">
			<date value="${taxon.date}" direction="forwards" units="years" uncertainty="0.0"/>
		</taxon>
</#list> 
	</taxa>
```

We can then replace the starting and data trees with the following code which will accept a tree from 
commandline and insert it into the xml.
```xml
    <rescaledTree id="startingTree" height="1.0">
        <newick usingDates="true">
            ${tree}
        </newick>
    </rescaledTree>
	
	<newick id="dataTree" usingDates="false" usingHeights="true">
		${tree}
	</newick>
```
And we can save it as SARS-CoV-2.mock.template

The xml can be generated with
```bash
beastgen -date_order -1 -date_prefix "|" -date_precision ./InprogressFiles/1K_SARS-CoV-2.mock.template  ./data/1K_SARS-CoV-2.nexus 1K_SARS-CoV-2_SG-thorney.xml 
```

This has the added benefit of reading the taxa out of the nexus file and by passes the need for a fasta file.


The template, `SARS-CoV-2.mock.template`, is small and easy to edit. It can also be used to generate similar analysis on much
larger datsets without the tedious task of editing a large xml. The operator weights are experimental and may need to be 
adjusted. 

A brief introduction to BEASTGen can be found [here](http://beast.community/beastgen). The Apache FreeMarker template
engine is quite powerful, which makes BEASTGen a useful tool for managing large datasets and analyses. 

## References
Zuckerkandl E , Pauling L. 1962. Molecular disease, evolution, and genic heterogeneity. In: Kasha M , Pullman B, editors. Horizons Biochem. New York: Academic Press. p. 189–222

[Xavier Didelot, Nicholas J Croucher, Stephen D Bentley, Simon R Harris, Daniel J Wilson, Bayesian inference of ancestral dates on bacterial phylogenetic trees, Nucleic Acids Research, Volume 46, Issue 22, 14 December 2018, Page e134](https://doi.org/10.1093/nar/gky783)

[Xavier Didelot, Igor Siveroni, Erik M Volz, Additive Uncorrelated Relaxed Clock Models for the Dating of Genomic Epidemiology Phylogenies, Molecular Biology and Evolution](https://doi.org/10.1093/molbev/msaa193)

## Help and documentation

The BEAST website: [http://beast.community](http://beast.community)

Tutorials: [http://beast.community/tutorials](http://beast.community/tutorials)

Frequently asked questions: [http://beast.community/faq](http://beast.community/faq)

