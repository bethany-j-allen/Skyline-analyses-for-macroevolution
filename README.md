---
author: Bethany J. Allen
level: Intermediate
title: Skyline analyses for macroevolution
subtitle: Applying BDSKY to species phylogenies
beastversion: 2.6.6
---


# Background

Bayesian phylodynamics uses the shape of a phylogenetic tree to infer characteristics of the population described by the phylogeny. Although widely applied to epidemiological datasets (see Skyline plots tutorial), the approach is yet to be used widely in macroevolution. In this case, skyline methods can be used to estimate parameters such as speciation, extinction and sampling rates over time, as well as the total number of lineages (usually species diversity). In this tutorial, we demonstrate how to apply the exponential coalescent and fossilised-birth-death skyline models, which both estimate piecewise-constant evolutionary rates through time, to a dinosaur supertree. The models differ in the temporal direction in which they are applied, and the assumptions they make about how the phylogeny is sampled.

----

# Programs used in this Exercise 

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

[BEAST2] (http://www.beast2.org) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014 --file Tutorial-Template/master-refs.bib %}. This tutorial uses the BEAST2 version 2.6.6.

### BEAUti - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs, and the interface will be the same, on all computing platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### Any programmer-friendly text editor

We will need to edit the XML files produced by BEAUti, for which we'll need a text editor. It's best to use one designed for programmers as these include nice features such as colouring code, which makes it more reader-friendly.

### Tracer

[Tracer] (http://beast.community/tracer) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v1.7.x.

### R / RStudio

We will be using R to analyze the output of the Birth-Death Skyline plot. RStudio provides a user-friendly graphical user interface to R that makes it easier to edit and run scripts. (It is not necessary to use RStudio for this tutorial.)

----

# Practical: Skyline analyses for macroevolution

In this tutorial we will estimate diversification rates for dinosaurs using a previously published supertree.

The aim of this tutorial is to:
- Learn how to set up a skyline analysis using a previously made phylogeny;
- Develop skills in simple XML hacking;
- Highlight the differences between exponential coalescent and fossilised-birth-death skylines.

## The data
We will be inferring our skyline parameters using a ready-made phylogeny containing 420 dinosaur species, published by {% cite Lloyd2008 --file Tutorial-Template/master-refs.bib %}. This phylogeny is a 'supertree', created using an informal method to collate several smaller dinosaur phylogenies into a larger one. Supertrees are typically cladograms, and must be timescaled using fossil data to estimate their branch lengths. The branch lengths we use here were inferred by {% cite Sakamoto2016 --file Tutorial-Template/master-refs.bib %}: they used fossil occurrences from the [Paleobiology Database] (http://paleobiodb.org) to infer the midpoint of the temporal range of each tip, and timescaled the phylogeny using the "equal" method, which distributes time evenly between the available branches.

## Setting up the Fossilised-Birth-Death skyline analysis
Many of the features we will need in our XML files are currently not implemented in BEAUti. However, we will start our analyses by creating XML files in BEAUti which will then serve as a template for us to alter by hand ("hack") later.

### Install BEAST2 packages
The coalescent-based skyline model is included in the core of BEAST2, but we need to install the **BDSKY** package, which contains the birth-death skyline model. We will also need the **feast** package, which will allow us to integrate some more complex features into our analyses. Installation of packages is done using the package manager, which is integrated into BEAUti.

>Open the **BEAST2 Package Manager** by navigating to **File > Manage Packages**.
>
>Install the **BDSKY** and **feast** packages by selecting them and clicking the **Install/Upgrade** button one at a time.

After the installation of a package, the program is on your computer, but BEAUti is unable to load the template files for the newly installed packages unless it is restarted.

>Close the **BEAST2 Package Manager** and **_restart_** BEAUti to fully load the **BDSKY** and **feast** packages.

### Creating the Analysis Files with BEAUti
The first step in setting up an analysis in BEAUti is to upload an alignment. In our case, we do not have an alignment as we are using a ready-made phylogeny instead, but we will still need to upload an alignment in order to initialise our XML. We have provided a dummy `.nexus` file for this purpose, which we will replace later with code that reads in our phylogeny.

>In the **Partitions** panel, import the nexus file with the empty alignment by navigating to **File > Import Alignment** in the menu and then finding the `empty.nexus` file on your computer, *or* drag and drop the file into the **BEAUti** window.

The empty alignment should appear in the **Partitions** tab, named `empty` and containing a single nucleotide site for a single taxon.

Most of the default tabs in BEAUti relate to inferring the phylogeny, which we will not be doing, so we can skip straight to the **Priors** tab. It's worth remembering that a lot of the default parameters in BEAUti are set up for analyses very different to ours, so it's worth checking these thoroughly.

>Select the **Priors** tab and choose **Birth Death Skyline BDSParam** as the tree prior.

The "standard" `contemporary` birth death skyline is parameterised using a reproductive number, a "become uninfectious" rate and a sampling proportion. However, here we are using the `BDSParam` model, which refers to parameterisation using a birth rate (here, speciation), a death rate (here, extinction) and an extant sampling proportion. As our phylogeny is of non-avian dinosaurs, for which we only have fossils and no extant (genetic) samples, we will later exchange the extant sampling proportion, usually denoted using {% eqinline rho %}, for a fossil sampling rate, usually denoted using {% eqinline psi %}.

Choosing sensible priors for these parameters is not straightforward, so we will select priors which are relatively unrestrictive. For the birth and death rates, we will use exponential priors with a mean of 1.0; this places more probability on small rates, but still permits rates which are towards the higher end of those estimated from living animals and plants {% cite HenaoDiaz2019 --file Tutorial-Template/master-refs.bib %}. For the sampling rate, we will also choose an exponential prior, this time with a mean of 0.2.

On the **initial =** buttons you will see two sets of square brackets. The first indicates what the starting value for that parameter will be, meaning its value in the first iteration of your chain. We need to alter our initialisation values to ensure that they sit within our prior distributions. The second contains two values which denote the limits of the range of values that our parameters are permitted to take, which are also important to consider carefully. Priors influence the probability of certain values being tested in your chain, but values which are improbable under your prior can still be selected if the signal in your data is strong enough; setting this range provides hard limits to your parameter values regardless of your priors. Our parameters are all rates, expressed per branch per million years, so the full set of values they can take ranges between 0.0 and infinity.

This is also the point where we express how many sections (here called **dimensions**) we want in our **piecewise constant** rates. Our rates will be assumed to be constant within these sections, but will be permitted to change at the break points between them. To keep our analysis simple (in the hope of a timely convergence!), we are going to give our rates **four** dimensions, corresponding to the major geological intervals spanned by the phylogeny: the Triassic, Jurassic, Early Cretaceous and Late Cretaceous.

>Change the **birth rate** prior from a uniform distribution to an exponential. Using the drop-down arrow on the left, check that the **mean** is set to 1.0. Click on the **initial =** button and change the **initialisation value** to 1.0. Check that the lower value is 0.0 and the upper value is `Infinity`. Change the **Dimension** to 4. Repeat these four steps for the **death rate** prior.
>
>Change the **rho** (sampling) prior from a uniform distribution to an exponential. Using the drop-down arrow on the left, change the **mean** value to 0.2. Click on the **initial =** button and change the **initialisation value** to 0.5. Check that the lower value is 0.0, and change the upper value to `Infinity`. Change the **Dimension** to 4.

We can leave the rest of the tabs as they are and save the XML file. We want to shorten the chain length and decrease the sampling frequency so the analysis completes in a reasonable time and the output files stay small. Note that we won't alter the **treelog** because we won't be using it, as our tree is fixed. (Keep in mind that it will be necessary to run a longer chain for parameters to mix properly and converge.)

>Navigate to the **MCMC** panel.
>
>Change the **Chain Length** from 10’000’000 to 1’000’000.
>
>Click on the arrow next to **tracelog** and change the **File Name** to `$(filebase).log` and set the Log Every to 1’000.
>
>Leave all other settings at their default values and save the file as `dinosaur_BDSKY.xml`.

We can now start "hacking" our XML template to remove the content we don't need and add some additional features.

>Open `dinosaur_BDSKY.xml` in your preferred text editor.

The first line sets out information about the format of the xml and its contents, which we can ignore. The next section is labelled `data`, and contains our dummy alignment. We don't need this section, so it can be commented out or simply deleted.

>Remove the `data` section of the XML:
>
>```xml
>    <data
> id="empty"
> spec="Alignment"
> name="alignment">
>                     <sequence id="seq_Tyrannosaurus_rex" spec="Sequence" taxon="Tyrannosaurus_rex" totalcount="4" value="N"/>
>                 </data>
>```

In it's place, we instead need to read in our phylogeny. We will do this by pasting in a new section, `tree`, which will use `TreeFromNewickFile` in the package **feast** to retrieve the phylogeny from our tree file, `Lloyd.tree`.

>Where the `data` section was in the XML, paste in the `tree` section:
>
>```xml
>  <tree id="tree"
>        spec="feast.fileio.TreeFromNewickFile" fileName="Lloyd.tree"
>        IsLabelledNewick="true" adjustTipHeights="false" />
>```

The section after this contains a block of statements labelled `map`, which links the different distributions we might use (for example, for our priors) to the BEAST2 code which defines them. We will leave this alone, and move on to the largest and final block, `run`, which describes the analysis we want to carry out.

The first subsection within the `run` block is labelled `state`. This describes the objects inferred within the analysis. The first object described is the `tree`, based on our dummy alignment. As we are not inferring our tree, it can be removed from this section.

>Remove the `tree` part of the `state` subsection of the XML:
>
>```xml
>        <tree id="Tree.t:empty" spec="beast.evolution.tree.Tree" name="stateNode">
>            <taxonset id="TaxonSet.empty" spec="TaxonSet">
>                <alignment idref="empty"/>
>            </taxonset>
>        </tree>
>```

Following this we see a description of our three inferred rates: birth, death and sampling. These lines include a lot of information, much of which we specified earlier via the **Priors** tab in BEAUti. As previously mentioned, our sampling parameter is currently set up as the extant sampling proportion, `rho`, which we want to exchange for a fossil sampling rate, `sampling`. For the sake of readability, we will change the `id` (name) of our sampling parameter here from `rhoBDS` to `samplingBDS`. Here, this is simply a label and does not change what the parameter does, but it is how the parameter is referenced in other parts of the XML, where its name will also need to be changed.

>Replace all instances of `rho` in the XML with `sampling`. There should be eight. This is most easily (and reliably) done using `Find & Replace` functionality.

To our three rates we will also add a fourth parameter, `origin`. This describes the time at which the most recent common ancestor of the clade diverged into its first two daughter species, i.e. the start of our evolutionary processes. Note that we want to infer this time because it is only the same as the timing of the first divergence in our sampled phylogeny if both of those initial species are represented in its tips. If this is not the case, the "true" origin must lie an unknown amount of time earlier than this first observed divergence.

In this section we need to set the limits of our origin time and its initial value. BEAST2 assumes that time runs from the present backwards (as is the case in our **Exponential Coalescent** model), so we will set our origin's lower limit to 0 (the age of the youngest tip), and to be maximally conservative, its upper limit to `Infinity`. We will set the starting value to 200 (200 million years before the youngest tip), which corresponds to an age of approximately 200 + 66 = 266Ma (if our youngest tip lies at the Cretaceous-Paleogene boundary, when non-avian dinosaurs are thought to have become extinct). This would place the origin of dinosaurs during the middle Permian, which is perhaps a little earlier than most palaeontologists would speculate (late Permian to Early Triassic, e.g. ADD CITATION), but is certainly adequate for our first iteration.

>Paste in the `origin` line to the end of the `parameter`block:
>
>`<parameter id="origin.t:empty" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="Infinity">200</parameter>`
>
>The `state` block should now look like this:
>
>```xml
><state id="state" spec="State" storeEvery="5000">
>        <parameter id="birthRateBDS.t:empty" spec="parameter.RealParameter" dimension="4" lower="0.0" name="stateNode" upper="Infinity">1.0</parameter>
>        <parameter id="deathRateBDS.t:empty" spec="parameter.RealParameter" dimension="4" lower="0.0" name="stateNode" upper="Infinity">1.0</parameter>
>        <parameter id="samplingBDS.t:empty" spec="parameter.RealParameter" dimension="4" lower="0.0" name="stateNode" upper="Infinity">0.5</parameter>
>	 <parameter id="origin.t:empty" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="Infinity">200</parameter>
></state>
>```

Note that each of our rate parameters have 4 **dimensions** but we do not specify this value for our origin; while our **piecewise constant** rates can change three times in each iteration (resulting in four fixed values), only a single origin value is needed, and we therefore do not need to specify its number of dimensions (as 1).

The next block, `init`, describes the initialisation processes, namely construction of the starting tree and population model, neither of which we need.

>Remove the `init` subsection of the XML:
>
>```xml
>    <init id="RandomTree.t:empty" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:empty" taxa="@empty">
>        <populationModel id="ConstantPopulation0.t:empty" spec="ConstantPopulation">
>            <parameter id="randomPopSize.t:empty" spec="parameter.RealParameter" name="popSize">1.0</parameter>
>        </populationModel>
>    </init>
>```

The next block of the XML is arguably the most important: it defines the **prior and posterior distributions** of our analysis, including the **model** we are using. The model itself lies in the third nested `distribution` line, labelled `BirthDeathSkyContemporaryBDSParam`. With the addition of line breaks to aid readability, it should currently look like this:

```xml
	<distribution id="BirthDeathSkyContemporaryBDSParam.t:empty"
                          spec="beast.evolution.speciation.BirthDeathSkylineModel"
                          birthRate="@birthRateBDS.t:empty"
                          conditionOnRoot="true"
                          contemp="true"
                          deathRate="@deathRateBDS.t:empty"
                          sampling="@samplingBDS.t:empty"
                          tree="@Tree.t:empty">
```

Although perhaps initially intimidating, we will walk through the various arguments here and what they mean. The model will also require some tweaking to fit our intended usage.

The first thing to note is that we're linking the parameters we defined earlier into the model. For example, `birthRate="@birthRateBDS.t:empty"` is simply telling the model the name we have given to our birth rate parameter. This is also true for `deathRate`. If you successfully changed all of the `rho` references to `sampling` earlier, you should also see that `sampling="@samplingBDS.t:empty"`. The name of the argument in the model is actually `samplingRate`, so we will correct that now.

>Change `sampling="@samplingBDS.t:empty"` to `samplingRate="@samplingBDS.t:empty"`.
		
We also need to link the origin parameter (which we defined earlier) into the model.
		
>Somewhere within the model line, add `origin="origin.t:empty"`.
		
The next argument to note is `conditionOnRoot="true"`. Whenever we conduct a Bayesian birth-death phylogenetic or phylodynamic analysis, we have to choose a single facet of the model which will remain fixed in some way, for all of the other parameters to move around. We describe the model as being **conditioned** on this value. Here we can see that the default value to condition on is the age of the origin (or root); as discussed earlier, this means assuming that the first divergence in our sampled phylogeny represents the "true" origin time. This might be a reasonable assumption for a well-sampled phylogeny which ideally incorporates both fossil and extant taxa, but that is not the case here, which is why we have instead chosen to infer our origin time. We must therefore choose something else to condition the model on.
		
Another option is to condition on ensuring that the model produces at least one sample (is observable). Depending on the type of sampling used, we can either use `conditionOnRhoSampling` or `conditionOnSurvival`: the former assumes at least one extant sample, while the latter assumes at least one sample of either type, be it an extant sample or an observed fossil. As our phylogeny of non-avian dinosaurs only includes extinct species, we will set these to "false" and "true" respectively.
		
>Change `conditionOnRoot="true"` to `conditionOnSurvival="true" conditionOnRhoSampling="false"`.
		
After this is the argument `contemp="true"`. This describes whether extant samples were collected in the present, i.e. at the age of the youngest tips. As we do not have extant samples, this is not relevant to our analysis, and can be removed. However, for our fossil sampling rate, we instead need to specify a different parameter, the `removalProbability`. When applying birth-death models to epidemiological datasets, sampling events are usually assumed to be associated with removal from the dataset: once a patient has been diagnosed with a pathogen, it is assumed that their pathogen is sampled and sequenced once (to enable inclusion in the phylogeny), and that they subsequently start receiving treatment, removing them from the population of diseased individuals. For our dataset, this would be equivalent to assuming that fossils are only deposited at the time of extinction. It would prevent multiple fossils existing for any single species, and mean than species for which fossils have been found could not be direct ancestors of other species in the phylogeny. While the relevance of these assumptions to supertrees is debatable, it is most logical to facilitate "sampled ancestors" in the case of fossil phylogenies (see {% cite Gavryushkina2014 --file Tutorial-Template/master-refs.bib %}), so we will set our `removalProbability` to 0.

>Change `contemp="true"` to `removalProbability="0.0"`.

The last argument currently in the model line links the tree into the model. We need to update the name of the tree object to correspond to our fixed phylogeny which we read in at the start of the XML.

>Change `tree="@Tree.t:empty"` to `tree="@tree"`.

The model is almost ready, but we want to add one last set of arguments to it. Our rates are **piecewise constant**, and we have already specified in the `state` block how many **dimensions** (constant intervals) each of our rates will have. By default, the break points in our rates are evenly spaced between the origin and the youngest tip, meaning that all of the intervals are the same length. However, we can add arguments to our model specifying when we would like our break points to be. For us, it makes sense to place these break points at the boundaries of geological intervals, so that we estimate a rate for each interval. As mentioned before, we are going to align our break points to the boundaries between the Triassic, Jurassic, Early Cretaceous and Late Cretaceous. To do this, we need to supply vectors of these times relative to the phylogeny. We are assuming that the youngest tip, at which `t=0`, lies at the Cretaceous-Paleogene boundary, so the vectors describe the cumulative duration of these geological intervals relative to this boundary. We need to specify a separate vector for each of our rates.

>Add to the end of the model line:
>
>```xml
>birthRateChangeTimes="0 32.55 77.05 133.35"
>deathRateChangeTimes="0 32.55 77.05 133.35"
>samplingRateChangeTimes="0 32.55 77.05 133.35"
>```
	
Beneath our model, you will see a line which sets the `samplingRate` to 0. This refers to fossil sampling, which we have now introduced into our model, so this line needs to be removed.

In its place, we instead need to add an instruction to the XML about the direction of time. BEAST2 assumes that time runs from the present backwards, as is the case in coalescent models. But in birth-death models, time instead runs forward, and we need to specify this using `reverseTimeArrays`. Here, we are telling the model to use the opposite direction of time to the conventional one in five dimensions: in order, this refers to our three rate parameters (`birth`, `death` and `sampling`) plus `rho` and the `removalProbability` (which we do not need but will specify anyway).

>Remove the line fixing `samplingRate` to 0:
>
>`<parameter id="samplingRateBDS.t:empty" spec="parameter.RealParameter" name="samplingRate">0.0</parameter>`
>
>Replace it with an instruction to `reverseTimeArrays`:
>
>`<reverseTimeArrays id="BooleanParameter.0" spec="parameter.BooleanParameter" dimension="5">true true true true true</reverseTimeArrays>`

And with that our model is complete! In summary, the model should now read (with arguments in any order):

```xml
	<distribution id="BirthDeathSkyContemporaryBDSParam.t:empty"
                          spec="beast.evolution.speciation.BirthDeathSkylineModel"
			  origin="origin.t:empty"
                          birthRate="@birthRateBDS.t:empty"
			  deathRate="@deathRateBDS.t:empty"
                          samplingRate="@samplingBDS.t:empty"
                          conditionOnSurvival="true"
			  conditionOnRhoSampling="false"
                          removalProbability="0.0"
                          tree="@tree"
			  birthRateChangeTimes="0 32.55 77.05 133.35"
			  deathRateChangeTimes="0 32.55 77.05 133.35"
			  samplingRateChangeTimes="0 32.55 77.05 133.35">
	<reverseTimeArrays id="BooleanParameter.0" spec="parameter.BooleanParameter" dimension="5">true true true true true</reverseTimeArrays>
	</distribution>
```

The next part of our XML specifies the shape of our prior distributions. As long as you changed all of the `rho` references to `sampling` previously, we don't need to modify these any further; we provided all of the necessary information in BEAUti.

The last part of the `distribution` block is labelled the `likelihood`, and determines how well the inferred tree fits the data. Once again, we are using a fixed phylogeny, and can remove this from our XML.

>Remove the `likelihood` block from the XML:
>
>```xml
> <distribution id="likelihood" spec="util.CompoundDistribution" useThreads="true">
>             <distribution id="treeLikelihood.empty" spec="ThreadedTreeLikelihood" data="@empty" tree="@Tree.t:empty">
>                 <siteModel id="SiteModel.s:empty" spec="SiteModel">
>                     <parameter id="mutationRate.s:empty" spec="parameter.RealParameter" estimate="false" name="mutationRate">1.0</parameter>
>                     <parameter id="gammaShape.s:empty" spec="parameter.RealParameter" estimate="false" name="shape">1.0</parameter>
>                     <parameter id="proportionInvariant.s:empty" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
>                     <substModel id="JC69.s:empty" spec="JukesCantor"/>
>                 </siteModel>
>                 <branchRateModel id="StrictClock.c:empty" spec="beast.evolution.branchratemodel.StrictClockModel">
>                     <parameter id="clockRate.c:empty" spec="parameter.RealParameter" estimate="false" name="clock.rate">1.0</parameter>
>                 </branchRateModel>
>             </distribution>
>         </distribution>
>```

The penultimate block of the XML describes our **operators**. These define the moves that are used to propose new parameter values in the next iteration of the MCMC, so are fundamentally important to how our chain explores parameter space. Some are relevant to tree construction, so we can remove these.
		
>Remove the tree operators:
>
>```xml
> <operator id="BDSKY_contemp_bds_treeScaler.t:empty" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:empty" weight="3.0"/>
> <operator id="BDSKY_contemp_bds_treeRootScaler.t:empty" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:empty" weight="3.0"/>
> <operator id="BDSKY_contemp_bds_UniformOperator.t:empty" spec="Uniform" tree="@Tree.t:empty" weight="30.0"/>
> <operator id="BDSKY_contemp_bds_SubtreeSlide.t:empty" spec="SubtreeSlide" tree="@Tree.t:empty" weight="15.0"/>
> <operator id="BDSKY_contemp_bds_narrow.t:empty" spec="Exchange" tree="@Tree.t:empty" weight="15.0"/>
> <operator id="BDSKY_contemp_bds_wide.t:empty" spec="Exchange" isNarrow="false" tree="@Tree.t:empty" weight="3.0"/>
> <operator id="BDSKY_contemp_bds_WilsonBalding.t:empty" spec="WilsonBalding" tree="@Tree.t:empty" weight="3.0"/>
>```
		
The next three operators are all labelled `ScaleOperator`, and when fired, scale the values of our `death`, `sampling` and `birth` rates respectively. At the end of the lines, you can see a `weight` specified. This defines how often these operators fire relative to each other. As all three of these rates are important to us, and we want to ensure that we have explored their parameter space equally, we will make their weights equal.
	
>Change the `weight` of the `deathRateScaler`, `samplingScaler`, and `birthRateScaler` to 10.0.
		
One thing we are currently missing is a `ScaleOperator` for our additional parameter, the `origin`, so we will add one. The timing of the origin is less important to us than our three rates, so we will give its operator a smaller weight.
		
>Add an `originScaler` to the `operator` block:
>
>`<operator id="BDSKY_contemp_bds_originScaler.t:empty" spec="ScaleOperator" parameter="@origin.t:empty" weight="3.0"/>`
		
These operators are currently set to scale each of the **dimensions** in our rate vectors individually, but we can also add operators which scale all of the dimensions in our vectors in concert. This could prove useful in quickly determining the typical size of our rates. This can be done by setting `scaleAll="true"` in additional `ScaleOperators`. 
		
>Add some `scaleAll` operators to the `operator` block:
>
>```xml
> <operator id="BDSKY_contemp_bds_deathRateScalerAll.t:empty" spec="ScaleOperator" parameter="@deathRateBDS.t:empty" scaleAll="true" weight="10.0"/>
> <operator id="BDSKY_contemp_bds_samplingScalerAll.t:empty" spec="ScaleOperator" parameter="@samplingBDS.t:empty" scaleAll="true" weight="10.0"/>
> <operator id="BDSKY_contemp_bds_birthRateScalerAll.t:empty" spec="ScaleOperator" parameter="@birthRateBDS.t:empty" scaleAll="true" weight="10.0"/>
>```
		
Beneath the `ScaleOperators` you can see one more operator, an `UpDownOperator`. It has `up` and `down` arguments, which are specified as our `birth` and `death` rates respectively. As currently defined, when this operator fires, it increases the birth rate and decreases the death rate. This operator can be very useful in cases where we expect two of our parameters to be correlated with one another, as speciation and extinction rates often are {% cite HenaoDiaz2019 --file Tutorial-Template/master-refs.bib %}. However, the current parameterisation places a negative directionality on this relationship, whereas we would actually expect higher speciation rates to be matched by higher extinction rates. Instead, we are going to raise our birth and death rates in concert, and instead decrease our sampling rate. We will also increase the weight of this operator to match that of our `scaleOperators`.
		
>Change the `UpDownOperator` parameterisation from:
>
>```xml
> <operator id="BDSKY_contemp_bds_updownBD.t:empty" spec="UpDownOperator" scaleFactor="0.75" weight="2.0">
>        <up idref="birthRateBDS.t:empty"/>
>        <down idref="deathRateBDS.t:empty"/>
>    </operator>
>```
>
> to:
>
>```xml
> <operator id="BDSKY_contemp_bds_updownBD.t:empty" spec="UpDownOperator" scaleFactor="0.75" weight="10.0">
>        <up idref="birthRateBDS.t:empty"/>
>        <up idref="deathRateBDS.t:empty"/>
>	 <down idref="samplingBDS.t:empty"/>
>    </operator>
>```
		
And finally we reach the last block of the XML, which determines the **logs** outputted by our BEAST2 analysis. We can see a `tracelog`, which records our parameters to a log file, a `screenlog`, which provides settings on what is printed to the screen during the analysis, and a `treelog`, which records the trees sampled during our analysis. Our `treelog` would simply save identical versions of our input tree, so to save on file space we will remove it. We will also remove the `OperatorSchedule` - this can be useful for fine-tuning operators but we do not need it here.
		
>Remove the `treelog` and `OperatorSchedule`:
>
>```xml
>   <logger id="treelog.t:empty" spec="Logger" fileName="$(tree).trees" logEvery="1000" mode="tree">
>     	<log id="TreeWithMetaDataLogger.t:empty" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:empty"/>
>   </logger>
>   <operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>
>```
		
With that, our XML is ready! It's finally time to run the analysis in BEAST2.

### Setting up the Exponential Coalescent skyline analysis




-------

# Tutorial style guide

## Text styling

This is how to write _italic text_.

This is how to write **bold text**.

This is how to write **_bold and italic text_**.

Do text superscripts like this 7^th, x^2y or  x^(2y + 3z).


## Lists

### Unnumbered lists

- Lorem ipsum dolor sit amet, consectetur adipiscing elit.
- Integer pharetra arcu ut nisl mollis ultricies.
	- Fusce nec tortor at enim cursus dictum.
	- Phasellus nec urna quis velit eleifend convallis sodales nec augue.
- In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
- Nam vitae turpis eu lacus imperdiet mollis id at augue.
- Sed sed turpis ac dolor mollis accumsan.


### Numbered lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	1. Fusce nec tortor at enim cursus dictum.
	2. Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.

### Mixed lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	* Fusce nec tortor at enim cursus dictum.
	* Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.


## Figures


<figure>
	<a id="fig:example1"></a>
	<img style="width:25%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 1: This figure is 25% of the page width.</figcaption>
</figure>


<figure>
	<a id="fig:example2"></a>
	<img style="width:10%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 2: This figure is only 10% of the page width.</figcaption>
</figure>



# Code

A bit of inline monospaced font can be made `like this`. Larger code blocks can be made by using the code environment:

Java:

```java
public class HelloWorld {

    public static void main(String[] args) {
        // Prints "Hello, World" to the terminal window.
        System.out.println("Hello, World");
    }

}
```

XML:

```xml
	<BirthDeathSkylineModel spec="BirthDeathSkylineModel" id="birthDeath" tree="@tree" contemp="true">
	      <parameter name="origin" id="origin" value ="100" lower="0."/>    
	      <parameter name="R0" id="R0" value="2" lower="0." dimension ="10"/>
	      <parameter name="becomeUninfectiousRate" id="becomeUninfectiousRate" value="1" lower="0." dimension ="10"/>
	      <parameter name="samplingProportion" id="samplingProportion" value="0."/>
	      <parameter name="rho" id="rho" value="1e-6" lower="0." upper="1."/>
	</BirthDeathSkylineModel>
```

R:

```R
	> myString <- "Hello, World!"
	> print (myString)
	[1] "Hello, World!"
```

# Equations

Inline equations: {% eqinline \dot{x} = \sigma(y-x) %}

Displayed equations: 
{% eq \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right) %}



## Instruction boxes

Use block-quotes for step-by-step instruction that the user should perform (this will produce a framed box on the website):

> The data we have is not the data we want, and the data we need is not the data we have.
> 
> We can input **any** formatted text in here:
>
> - Even
> - Lists
>
> or equations:
>
> {% eq (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right) %}






# Hyperlinks

Add links to figures like this: 

- [Figure 1](#fig:example1) is 25% of the page width.
- [Figure 2](#fig:example2) is 10% of the page width. 

Add links to external URLs like [this](http://www.google.com). 

Links to equations or different sections within the same document are a little buggy.


----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Tutorial-Template/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Tutorial-Template/master-refs.bib %}

