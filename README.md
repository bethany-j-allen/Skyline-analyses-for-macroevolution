---
author: Bethany J. Allen
level: Intermediate
title: Skyline analyses for macroevolution
subtitle: Applying BDSKY to species phylogenies
beastversion: 2.6.6
---


# Background

Bayesian phylodynamics uses the shape of a phylogenetic tree to infer characteristics of the population represented in the phylogeny. Although widely applied to epidemiological datasets, the approach is yet to be used widely in macroevolution. In this case, skyline methods can be used to estimate parameters such as speciation, extinction and sampling rates over time, as well as the total number of lineages (usually species diversity). In this tutorial, we demonstrate how to apply the exponential coalescent and fossilised-birth-death skyline models, which both estimate piecewise-constant evolutionary rates through time, to a dinosaur supertree. The models differ in the temporal direction in which they are applied, and the assumptions they make about how the phylogeny is sampled. We recommend reading the **Skyline plots** tutorial before attempting this one, as it covers a lot of the theory behind the models which we will not repeat here, except to highlight points which are relevant to macroevolutionary datasets.

----

# Programs used in this Exercise 

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

[BEAST2](http://www.beast2.org) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014 --file Tutorial-Template/master-refs.bib %}. This tutorial uses the BEAST2 version 2.6.6.

### BEAUti - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs, and the interface will be the same, on all computing platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### Any programmer-friendly text editor

We will need to edit the XML files produced by BEAUti, for which we'll need a text editor. It's best to use one designed for programmers as these include nice features such as colouring code, which makes it more reader-friendly.

### Tracer

[Tracer](http://beast.community/tracer) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v1.7.x.

### R / RStudio

We will be using R to analyze the output of the Birth-Death Skyline plot. RStudio provides a user-friendly graphical user interface to R that makes it easier to edit and run scripts. (It is not necessary to use RStudio for this tutorial.)

----

# Practical: Skyline analyses for macroevolution

In this tutorial we will estimate diversification rates for dinosaurs using a previously published supertree. The evolutionary history of dinosaurs is somewhat controversial; although non-avian dinosaurs are well agreed to have become extinct as a result of an asteroid impact at the Cretaceous-Paleogene boundary, it has been fiercely debated whether the clade was already in decline prior to this event (e.g. {% cite Brusatte2015 --file Tutorial-Template/master-refs.bib %}, {% cite Benson2018 --file Tutorial-Template/master-refs.bib %}).

The aim of this tutorial is to:
- Learn how to set up a skyline analysis using a previously made phylogeny;
- Develop skills in XML hacking;
- Highlight the differences between exponential coalescent and fossilised-birth-death skylines.

If this tutorial has been useful to you, please consider citing {% cite Allen2023 --file Tutorial-Template/master-refs.bib %}, from which these analyses are taken.

## The data
We will be inferring our skyline parameters using a ready-made phylogeny containing 420 dinosaur species, published by {% cite Lloyd2008 --file Tutorial-Template/master-refs.bib %}. This phylogeny is a 'supertree', created using an informal method to collate several smaller dinosaur phylogenies into a larger one. Supertrees are typically cladograms, and must be timescaled using fossil data to estimate their branch lengths. The branch lengths we use here were inferred by {% cite Sakamoto2016 --file Tutorial-Template/master-refs.bib %}: they used fossil occurrences from the [Paleobiology Database](http://paleobiodb.org) to infer the midpoint of the temporal range of each tip, and timescaled the phylogeny using the "equal" method, which distributes time evenly between the available branches.

## Install BEAST2 packages
The coalescent-based skyline model is included in the core of BEAST2, but we need to install the **BDSKY** package, which contains the birth-death skyline model. We will also need the **feast** package, which will allow us to integrate some more complex features into our analyses. Installation of packages is done using the package manager, which is integrated into BEAUti.

>Open the **BEAST2 Package Manager** by navigating to **File > Manage Packages**.
>
>Install the **BDSKY** and **feast** packages by selecting them and clicking the **Install/Upgrade** button one at a time.

After the installation of a package, the program is on your computer, but BEAUti is unable to load the template files for the newly installed packages unless it is restarted.

>Close the **BEAST2 Package Manager** and **_restart_** BEAUti to fully load the **BDSKY** and **feast** packages.

## Setting up the Exponential Coalescent Skyline analysis
We will start with the simpler of the two models, the **exponential coalescent**. This model is similar in construction to the **constant coalescent** model applied in the **Skyline plots** tutorial, but instead of assuming that the size of the population remains constant during our individual time intervals, we instead assume that they are experiencing **exponential growth or decline**. The advantage of using this model is that while the constant coalescent estimates **constant effective population sizes** within each time interval, the exponential coalescent instead estimates **diversification rates** for each of the time intervals, which is what we would like to infer for our dinosaurs.

Many of the features we will need in our XML files are not currently implemented in BEAUti. However, for both models, we will start our analyses by creating XML files in BEAUti which will then serve as a template for us to alter by hand ("hack") later. 

### Creating the Analysis Files with BEAUti
The first step in setting up an analysis in BEAUti is to upload an alignment. In our case, we do not have an alignment as we are using a ready-made phylogeny instead, but we will still need to upload an alignment in order to initialise our XML. We have provided a dummy `.nexus` file for this purpose, which we will replace later with code that reads in our phylogeny.

>In the **Partitions** panel, import the nexus file with the empty alignment by navigating to **File > Import Alignment** in the menu and then finding the `empty.nexus` file on your computer, *or* drag and drop the file into the **BEAUti** window.

The empty alignment should appear in the **Partitions** tab, named `empty` and containing a single nucleotide site for a single taxon.

Most of the default tabs in BEAUti relate to inferring the phylogeny, which we will not be doing, so we can skip straight to the **Priors** tab. It's worth remembering that a lot of the default parameters in BEAUti are set up for analyses very different to ours, so it's worth checking these thoroughly.

>Select the **Priors** tab and choose **Coalescent Exponential Population** as the tree prior.

Here we see that the model has two parameters, `ePopSize` and `growthRate`. `ePopSize` refers to the **effective population size** at the start of the coalescent process. Because coalescent models consider time from the present backwards (see **Skyline plots** tutorial), this therefore refers to the size of the population at the end of our youngest time interval. The tips in our phylogeny are species, and so in this context, our effective population size can be considered to be analogous to **total species richness**. Our `growthRate` is simply our **diversification rate**.

We actually only need the population size prior, and will remove the growth rate prior later. The default `ePopSize` prior is a 1/X (or `OneOnX`), which is a good choice when you have little knowledge about what the shape of your prior should be, placing reduced probability on higher values. We will keep the shape of this prior as the default, but change the initialisation values (the value of the parameters in their first iteration of the chain). We will set the starting `ePopSize` to 1, and the `growthRate` to 0.

>Click on the **initial =** button for `ePopSize` and change the **initialisation value** to 1.0.
>
>Click on the **initial =** button for `growthRate` and change the **initialisation value** to 0.0.

We will go into much more detail on the contents of this tab later, when setting up our birth-death model.

We can leave the rest of the tabs as they are and save the XML file. We want to shorten the chain length and decrease the sampling frequency so the analysis completes in a reasonable time and the output files stay small. Note that we won't alter the **treelog** because we won't be using it, as our tree is fixed. (Keep in mind that it will be necessary to run a longer chain for parameters to mix properly and converge.)

>Navigate to the **MCMC** panel.
>
>Change the **Chain Length** from 10’000’000 to 1’000’000.
>
>Click on the arrow next to **tracelog** and change the **File Name** to `$(filebase).log` and set the Log Every to 1’000.
>
>Leave all other settings at their default values and save the file as `dinosaur_coalescent.xml`.

### Amending the Analysis Files
We can now start "hacking" our XML template to ensure that it includes our desired features, including some that are not available in BEAUti.

>Open `dinosaur_coal.xml` in your preferred text editor.

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

After this we see descriptions of our two parameters, `ePopSize` and `growthRate`. Rather than one growth rate, we actually want to infer one for each of our **time intervals** of interest. To keep our analysis simple (in the hope of a timely convergence!), we are going to use **four** time bins, corresponding to the major geological intervals spanned by the phylogeny: the Triassic, Jurassic, Early Cretaceous and Late Cretaceous. We will therefore copy and paste the `growthRate` parameter three times, giving each a new name (`id`).

>Convert the single `growthRate` parameter into four with different `ids`:
>
>```xml
> <state id="state" spec="State" storeEvery="5000">
> <parameter id="ePopSize.t:empty" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
> <parameter id="growthRate1" spec="parameter.RealParameter" name="stateNode">0.0</parameter>
> <parameter id="growthRate2" spec="parameter.RealParameter" name="stateNode">0.0</parameter>
> <parameter id="growthRate3" spec="parameter.RealParameter" name="stateNode">0.0</parameter>
> <parameter id="growthRate4" spec="parameter.RealParameter" name="stateNode">0.0</parameter>
> </state>
>```

The next block, `init`, describes the initialisation processes, namely construction of the starting tree and assoicated population model, neither of which we need.

>Remove the `init` subsection of the XML:
>
>```xml
>    <init id="RandomTree.t:empty" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:empty" taxa="@empty">
>        <populationModel id="ConstantPopulation0.t:empty" spec="ConstantPopulation">
>            <parameter id="randomPopSize.t:empty" spec="parameter.RealParameter" name="popSize">1.0</parameter>
>        </populationModel>
>    </init>
>```

The next block of the XML is arguably the most important: it defines the **prior and posterior distributions** of our analysis, including the **model** we are using. The model itself is described in the third nested `distribution` line, labelled `CoalescentExponential`. With the addition of line breaks to aid readability, it should currently look like this:

```xml
<distribution id="CoalescentExponential.t:empty" spec="Coalescent">
                <populationModel id="ExponentialGrowth.t:empty"
		spec="ExponentialGrowth"
		growthRate="@growthRate.t:empty"
		popSize="@ePopSize.t:empty"/>
                <treeIntervals id="TreeIntervals.t:empty" spec="TreeIntervals" tree="@Tree.t:empty"/>
</distribution>
```

The arguments determine the model to be used, and link the two parameters in the model, `growthRate` and `popSize`, to our definitions of them elsewhere in the XML. We can see that at present, the model uses a single growth rate for the full duration of the coalescent process. Instead of using this model, we are going to use a variation on it which is available in the package **feast**. This model is described as a `CompoundPopulationModel`, and allows multiple intervals to be defined over the course of our coalescent process, each of which can have their own growth rate. This will enable us to estimate a different diversification rate for each of our four time intervals of interest. The basic parameterisation of the model, as described on the [feast Github page](https://github.com/tgvaughan/feast), is as follows:

```xml
<populationModel spec="CompoundPopulationModel">
    <populationModel spec="ConstantPopulation"> <popSize spec="RealParameter" value="5.0"/></populationModel>
    <populationModel spec="ConstantPopulation"> <popSize spec="RealParameter" value="10.0"/></populationModel>
    <populationModel spec="ConstantPopulation"> <popSize spec="RealParameter" value="2.0"/></populationModel>
    <changeTimes spec="RealParameter" value="1.0 3.0"/>
</populationModel>
```

We will be changing the `ConstantPopulation` specification to `ExponentialGrowth`, to match our current model. We will also link the size of the population at the end of one interval as the starting population size in the next, using the `makeContinuous="true"` argument.	
	
As well as defining models for three time intervals here, you can also see an object called `changeTimes`. Here we can provide a vector which describes when the boundaries between our time intervals should be, so it is easy for us to set these as our geological interval boundaries. The vector needs to provide the boundaries relative to the timescale of the phylogeny. We are assuming that the youngest tip lies at the Cretaceous-Paleogene boundary, so the times will need to describe the cumulative duration of these geological intervals relative to this boundary.

>Replace the current `populationModel` with the following:
>
>```xml
> <populationModel spec="feast.popmodels.CompoundPopulationModel" makeContinuous="true">
>           <populationModel spec="ExponentialGrowth" popSize="@ePopSize.t:empty" growthRate="@growthRate1"/>
>           <populationModel spec="ExponentialGrowth" popSize="1.0" growthRate="@growthRate2"/>
>           <populationModel spec="ExponentialGrowth" popSize="1.0" growthRate="@growthRate3"/>
>           <populationModel spec="ExponentialGrowth" popSize="1.0" growthRate="@growthRate4"/>
>           <changeTimes spec="parameter.RealParameter" value="32.55 77.05 133.35"/>
> </populationModel>
>```
	
We also need to change the `tree` name in `treeIntervals` so that it links with the phylogeny we are reading in at the start of the XML.
	
>Change `tree="@Tree.t:empty"` to `tree="@tree"` in the `treeIntervals` line.
	
As the fit of the model to the phylogeny is what we are interested in, we are going to rename this small block `likelihood` to ensure that this is what is recorded in our output logs.
	
>Change `id="CoalescentExponential.t:empty"` to `id="likelihood"`.
	
The model block should now look like this:
	
```xml
<distribution id="likelihood" spec="Coalescent">
	<populationModel spec="feast.popmodels.CompoundPopulationModel" makeContinuous="true">
           <populationModel spec="ExponentialGrowth" popSize="@ePopSize.t:empty" growthRate="@growthRate1"/>
           <populationModel spec="ExponentialGrowth" popSize="1.0" growthRate="@growthRate2"/>
           <populationModel spec="ExponentialGrowth" popSize="1.0" growthRate="@growthRate3"/>
           <populationModel spec="ExponentialGrowth" popSize="1.0" growthRate="@growthRate4"/>
           <changeTimes spec="parameter.RealParameter" value="32.55 77.05 133.35"/>
	</populationModel>
        <treeIntervals id="TreeIntervals.t:empty" spec="TreeIntervals" tree="@tree"/>
</distribution>
```
	
In the next few lines of our XML we can see that the shape of our prior distributions are defined. As mentioned before, we will keep the prior for `ePopSize` as the default 1/X. However, we are defining the shape of our `growthRate` parameters within our model so we can remove this prior from the XML.

>To tidy up, cut the `<distribution id="prior" spec="util.CompoundDistribution">` line from above the model block and paste it below this block (immediately above the `ePopSize` prior).
>
>Remove the `growthRate` prior:
>
>```xml
> <prior id="GrowthRatePrior.t:empty" name="distribution" x="@growthRate.t:empty">
>                <OneOnX id="OneOnX.3" name="distr"/>
> </prior>
>```

The last part of the `distribution` block is labelled the `likelihood`, and contains the `treeLikelihood`, which determines how well the inferred tree fits the data. As we are using a fixed phylogeny, we can remove this from our XML (we have already renamed our model the `likelihood` anyway).

>Remove the `likelihood` block containing the `treeLikelihood`:
>
>```xml
> <distribution id="likelihood" spec="util.CompoundDistribution" useThreads="true">
>   <distribution id="treeLikelihood.empty" spec="ThreadedTreeLikelihood" data="@empty" tree="@Tree.t:empty">
>     <siteModel id="SiteModel.s:empty" spec="SiteModel">
>       <parameter id="mutationRate.s:empty" spec="parameter.RealParameter" estimate="false" name="mutationRate">1.0</parameter>
>       <parameter id="gammaShape.s:empty" spec="parameter.RealParameter" estimate="false" name="shape">1.0</parameter>
>       <parameter id="proportionInvariant.s:empty" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
>       <substModel id="JC69.s:empty" spec="JukesCantor"/>
>     </siteModel>
>     <branchRateModel id="StrictClock.c:empty" spec="beast.evolution.branchratemodel.StrictClockModel">
>       <parameter id="clockRate.c:empty" spec="parameter.RealParameter" estimate="false" name="clock.rate">1.0</parameter>
>     </branchRateModel>
>   </distribution>
> </distribution>
>```

The penultimate block of the XML describes our **operators**. These define the moves that are used to propose new parameter values in the next iteration of the MCMC, so are fundamentally important to how our chain explores parameter space. Some are relevant to tree construction, so we can remove these.
		
>Remove the tree operators:
>
>```xml
> <operator id="CoalescentExponentialTreeScaler.t:empty" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:empty" weight="3.0"/>
> <operator id="CoalescentExponentialTreeRootScaler.t:empty" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:empty" weight="3.0"/>
> <operator id="CoalescentExponentialUniformOperator.t:empty" spec="Uniform" tree="@Tree.t:empty" weight="30.0"/>
> <operator id="CoalescentExponentialSubtreeSlide.t:empty" spec="SubtreeSlide" tree="@Tree.t:empty" weight="15.0"/>
> <operator id="CoalescentExponentialNarrow.t:empty" spec="Exchange" tree="@Tree.t:empty" weight="15.0"/>
> <operator id="CoalescentExponentialWide.t:empty" spec="Exchange" isNarrow="false" tree="@Tree.t:empty" weight="3.0"/>
> <operator id="CoalescentExponentialWilsonBalding.t:empty" spec="WilsonBalding" tree="@Tree.t:empty" weight="3.0"/>
>```

The next operator is a `ScaleOperator` which, when fired, **scales** the value of our `ePopSize` parameter. This is followed by a `RealRandomWalkOperator`, which alters our `growthRate` parameter via a **random walk** process. The first thing we will do is ensure that we have a `RealRandomWalkOperator` for each of our four growth rates.

>Copy and paste the `RealRandomWalkOperator` three times, and change the `id` and `parameter` names to our four separate growth rates, to create this:
>
>```xml
> <operator id="GrowthRateRandomWalk1.t:empty" spec="RealRandomWalkOperator" parameter="@growthRate1" weight="3.0" windowSize="1.0"/>
> <operator id="GrowthRateRandomWalk2.t:empty" spec="RealRandomWalkOperator" parameter="@growthRate2" weight="3.0" windowSize="1.0"/>
> <operator id="GrowthRateRandomWalk3.t:empty" spec="RealRandomWalkOperator" parameter="@growthRate3" weight="3.0" windowSize="1.0"/>
> <operator id="GrowthRateRandomWalk4.t:empty" spec="RealRandomWalkOperator" parameter="@growthRate4" weight="3.0" windowSize="1.0"/>
>```

One of the parameters in the `RealRandomWalkOperator` is the `windowSize`, which describes the permitted amplitude of parameter change. We will set this to 0.05 for each of our four rates.
	
>Change the `windowSize` of each of the `growthRate` operators to **0.05**.

You can also see that each of our operators has a `weight` specified. This defines how often these operators fire relative to each other. All five of the weights have a value of 3.0, meaning that they are equally likely to fire. This is a reasonable set-up for our analysis so we will leave these values as they are.

And finally we reach the last block of the XML, which determines the **logs** outputted by our BEAST2 analysis. We can see a `tracelog`, which records our parameters to a log file, a `screenlog`, which provides settings on what is printed to the screen during the analysis, and a `treelog`, which records the trees sampled during our analysis. Our `treelog` would simply save identical versions of our input tree, so to save on file space we will remove it. We will also remove the `OperatorSchedule` - this can be useful for fine-tuning operators but we do not need it here.

>Remove the `treelog` and `operatorschedule`:
>
>```xml
> <logger id="treelog.t:empty" spec="Logger" fileName="$(tree).trees" logEvery="1000" mode="tree">
>        <log id="TreeWithMetaDataLogger.t:empty" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:empty"/>
> </logger>
> <operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>
>```

We also need to remove parameters from the `tracelog` and `screenlog` that no longer exist in our XML, and ensure that we are logging all four of our `growthRate` parameters.
		
>Remove the `treeLikelihood`, `TreeHeight`, and `CoalescentExponential` parameters from the `tracelog`. Copy and paste <log idref="growthRate.t:empty"/> three times, changing the `idref` to `growthRate1`, `growthRate2`, `growthRate3`, and `growthRate4`.
		
With that, our XML is ready!

>**Save** all changes to your XML file and close it.

We can now run the analysis in BEAST2. It's important to have `Lloyd.tree` saved somewhere that BEAST2 can access it.

>Download and save `Lloyd.tree` in the same folder as your BEAST2 program. Open the program and select `dinosaur_coal.xml` (or `dinosaur_coal_final.xml` if you're using our ready-made version). If you have **BEAGLE** installed tick the box to **Use BEAGLE library if available**, which will make the run faster. Hit **Run** to start the analysis.
>
>**OR**
>
>Download and save `Lloyd.tree` in the folder containing your XML file. Find the BEAST2 executable in **BEAST_2.X.X** (depending on your version) **> bin**. Right-click on the **beast** executable and select **Create shortcut** on Windows or **Make alias** on Mac. Cut and paste the created shortcut/alias into the folder containing your analysis files. If you open your **terminal** and navigate to the folder containing your files, you should now be able to run the analysis through the terminal using `beast dinosaur_coal.xml` (or `beast dinosaur_coal_final.xml` if you're using our ready-made version).

The analysis should take about XX minutes to run. In the meantime, you can start setting up the fossilised-birth-death XML.

## Setting up the Fossilised-Birth-Death Skyline analysis
As with the exponential coalscent model, many of the features we will need in our XML file are not yet implemented in BEAUti, but we will start our analyses by creating XML files in BEAUti.

### Creating the Analysis Files with BEAUti
As before, we need to start by uploading our dummy `.nexus` file as an alignment.

>In the **Partitions** panel, import the nexus file with the empty alignment by navigating to **File > Import Alignment** in the menu and then finding the `empty.nexus` file on your computer, *or* drag and drop the file into the **BEAUti** window.

Once again, we will skip straight to the **Priors** tab.

>Select the **Priors** tab and choose **Birth Death Skyline BDSParam** as the tree prior.

The "standard" `contemporary` birth death skyline is parameterised using a reproductive number, a "become uninfectious" rate and a sampling proportion. However, here we are using the `BDSParam` model, which refers to parameterisation using a birth rate (here, speciation), a death rate (here, extinction) and an extant sampling proportion. As our phylogeny is of non-avian dinosaurs, for which we only have fossils and no extant (genetic) samples, we will later exchange the extant sampling proportion, usually denoted using {% eqinline rho %}, for a fossil sampling rate, usually denoted using {% eqinline psi %}.

Choosing sensible priors for these parameters is not straightforward, so we will select priors which are relatively unrestrictive. For the birth and death rates, we will use **exponential** priors with a mean of 1.0; this places more probability on small rates, but still permits rates which are towards the higher end of those estimated from living animals and plants {% cite HenaoDiaz2019 --file Tutorial-Template/master-refs.bib %}. For the sampling rate, we will also choose an exponential prior, this time with a mean of 0.2.

On the **initial =** buttons you will see two sets of square brackets. The first indicates what the starting value for that parameter will be; we need to alter our initialisation values to ensure that they sit within our prior distributions. The second contains two values which denote the limits of the range of values that our parameters are permitted to take, which are also important to consider carefully. Priors influence the probability of certain values being tested in your chain, but values which are improbable under your prior can still be selected if the signal in your data is strong enough. Setting this range provides hard limits to your parameter values regardless of your priors. Our parameters are all rates, expressed per branch per million years, so the full set of values they can take ranges between 0.0 and infinity.

This is also the point where we express how many sections (here called **dimensions**) we want in our **piecewise constant** rates. Our rates will be assumed to be constant within these sections, but will be permitted to change at the break points between them. Similar to the exponential coalescent model, we are going to give our rates **four** dimensions, corresponding to the Triassic, Jurassic, Early Cretaceous and Late Cretaceous.

>Change the **birth rate** prior from a uniform distribution to an exponential. Using the drop-down arrow on the left, check that the **mean** is set to 1.0. Click on the **initial =** button and change the **initialisation value** to 1.0. Check that the lower value is 0.0 and the upper value is `Infinity`. Change the **Dimension** to 4. Repeat these four steps for the **death rate** prior.
>
>Change the **rho** (sampling) prior from a uniform distribution to an exponential. Using the drop-down arrow on the left, change the **mean** value to 0.2. Click on the **initial =** button and change the **initialisation value** to 0.5. Check that the lower value is 0.0, and change the upper value to `Infinity`. Change the **Dimension** to 4.

We can leave the rest of the tabs as they are and save the XML file. We will again shorten the chain length and decrease the sampling frequency of our analysis, so the analysis completes in a reasonable time and the output files stay small.

>Navigate to the **MCMC** panel.
>
>Change the **Chain Length** from 10’000’000 to 1’000’000.
>
>Click on the arrow next to **tracelog** and change the **File Name** to `$(filebase).log` and set the Log Every to 1’000.
>
>Leave all other settings at their default values and save the file as `dinosaur_BDSKY.xml`.

### Amending the Analysis Files
It is now time to "hack" our XML template, to remove the content we don't need and add some additional features. The first few steps are the same as for our exponential coalescent XML.

>Open `dinosaur_BDSKY.xml` in your preferred text editor.
>
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
>
>Where the `data` section was in the XML, paste in the `tree` section:
>
>```xml
>  <tree id="tree"
>        spec="feast.fileio.TreeFromNewickFile" fileName="Lloyd.tree"
>        IsLabelledNewick="true" adjustTipHeights="false" />
>```
>
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
>`<parameter id="origin.t:empty" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="Infinity">200.0</parameter>`
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

As before, the `init` block is not needed for this analysis.

>Remove the `init` subsection of the XML:
>
>```xml
>    <init id="RandomTree.t:empty" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:empty" taxa="@empty">
>        <populationModel id="ConstantPopulation0.t:empty" spec="ConstantPopulation">
>            <parameter id="randomPopSize.t:empty" spec="parameter.RealParameter" name="popSize">1.0</parameter>
>        </populationModel>
>    </init>
>```

The next block describes the model, which is labelled `BirthDeathSkyContemporaryBDSParam`. With the addition of line breaks to aid readability, it should currently look like this:

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
		
>Somewhere within the model line, add `origin="@origin.t:empty"`.
		
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

In it's place, we instead need to add an instruction to the XML about the direction of time. BEAST2 assumes that time runs from the present backwards, as is the case in coalescent models. But in birth-death models, time instead runs forward, and we need to specify this using `reverseTimeArrays`. Here, we are telling the model to use the opposite direction of time to the conventional one in five dimensions: in order, this refers to our three rate parameters (`birth`, `death` and `sampling`) plus `rho` and the `removalProbability` (which we do not need but will specify anyway).

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
			  origin="@origin.t:empty"
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

The next part of our XML specifies the shape of our prior distributions. As long as you changed all of the `rho` references to `sampling` previously, we don't need to modify the priors we already have any further; we provided all of the necessary information in BEAUti. We do, however, need to add a prior on our `origin`. You may remember that when specifying the parameter, we set our starting value to 200 (corresponding to 266Ma, which is on the early side of expert estimates for the origin of dinosaurs). To keep things simple, we will set the origin prior to a **uniform** distribution between 0 and 200.

>Add an `origin` prior to the `priors` block:
>
>```xml
> <prior id="originPrior" name="distribution" x="@origin.t:empty">
>          <Uniform name="distr" upper="200"/>
> </prior>
>```

As before, the last part of the `distribution` block is labelled the `likelihood`, and contains the `treeLikelihood`, which we don't need because we are using a fixed phylogeny. What we consider the `likelihood` in our analysis is instead how well the model parameters fit our phylogeny, so this time we will cut and paste our model into this part of the XML.

>Remove the `treeLikelihood` block from the `likelihood` part of the XML:
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
>
>Cut and paste the model into this section and remove `useThreads="true"`:
>```xml
> <distribution id="likelihood" spec="util.CompoundDistribution">
>            <distribution id="BirthDeathSkyContemporaryBDSParam.t:empty"
>                            spec="beast.evolution.speciation.BirthDeathSkylineModel"
>                            origin="@origin.t:empty"
>                            birthRate="@birthRateBDS.t:empty"
>                            deathRate="@deathRateBDS.t:empty"
>                            samplingRate="@samplingBDS.t:empty"
>                            conditionOnSurvival="true"
>                            conditionOnRhoSampling="false"
>                            removalProbability="0.0"
>                            tree="@tree"
>                            birthRateChangeTimes="0 32.55 77.05 133.35"
>                            deathRateChangeTimes="0 32.55 77.05 133.35"
>                            samplingRateChangeTimes="0 32.55 77.05 133.35">
>                            <reverseTimeArrays id="BooleanParameter.0"
>                              spec="parameter.BooleanParameter"
>                              dimension="5">true true true true true</reverseTimeArrays>
>            </distribution>
> </distribution>

We now see the **operators**. As for the exponential coalescent model, we will remove the operators related to tree construction, which we don't need.
		
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
		
The next three operators are all labelled `ScaleOperator`, and when fired, scale the values of our `death`, `sampling` and `birth` rates respectively. As before, you can see a `weight` specified, defining how often these operators fire relative to each other. As all three of these rates are important to us, and we want to ensure that we have explored their parameter space equally, we will make their weights equal.
	
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
		
Once again, the last block of the XML determines the **logs** outputted by the analysis. We will remove the `treelog` and `operatorschedule` to save file space.
		
>Remove the `treelog` and `operatorschedule`:
>
>```xml
>   <logger id="treelog.t:empty" spec="Logger" fileName="$(tree).trees" logEvery="1000" mode="tree">
>     	<log id="TreeWithMetaDataLogger.t:empty" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:empty"/>
>   </logger>
>   <operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>
>```
		
We also need to remove parameters from the `tracelog` and `screenlog` that no longer exist in our XML, and log the `origin`.
		
>Remove the `treeLikelihood` and `TreeHeight` parameters from the `tracelog`, and add `<log idref="origin.t:empty"/>`. 
		
Our XML is finally ready!

>**Save** all changes to your XML file and close it.

We can now run the analysis in BEAST2. As with the exponential coalescent model, it's important to have `Lloyd.tree` saved somewhere that BEAST2 can access it (check back to the instructions for the running the **exponential coalescent** XML if you skipped this).

>Open the program and select `dinosaur_BDSKY.xml` (or `dinosaur_BDSKY_final.xml` if you're using our ready-made version). If you have **BEAGLE** installed tick the box to **Use BEAGLE library if available**, which will make the run faster. Hit **Run** to start the analysis.
>
>**OR**
>
>Ensure your XML file and `Lloyd.tree` are saved in the same folder, which also contains the BEAST2 shortcut/alias you created earlier. Open your **terminal** and navigate to the folder containing your files, then run the analysis through the terminal using `beast dinosaur_BDSKY.xml` (or `beast dinosaur_BDSKY_final.xml` if you're using our ready-made version).

The analysis should take about XX minutes to run.

##Visualising the results
Once the BEAST2 analyses have finished running, we will use **R** to plot our skylines. The log files are relatively easy to handle in R, so we have provided custom code for this rather than using an R package, although this code does require the **tidyverse** to be installed (specifically, we will use **dplyr** and **ggplot2**).

The first step is to install the `tidyverse` if you haven't previously, and then to load the package.

```R
install.packages("tidyverse")

library(tidyverse)
```

###The Exponential Coalescent model results
First we will take a look at our **exponential coalescent** skyline. The log file from this analysis needs to be read into R, either with or without setting the working directory (the location where R will look for your files). We will also immediately trim the log to remove the first 10% of iterations as burn-in.

```R
# Navigate to Session > Set Working Directory > Choose Directory (on RStudio)
# or change file name to the full path to the log file
#(Use "dinosaur_coal_final.log" if you used our pre-cooked XML)
coal_file <- "dinosaur_coal.log"

#Read in coalescent log and trim 10% burn-in
coalescent <- read.table(coal_file, header = T) %>% slice_tail(prop = 0.9)
```

The next job is to summarise the values across the remaining iterations in the log. We will pivot the table so that our diversification rate estimates sit in a single column, then estimate our median and 95% **highest posterior density** (HPD; Bayesian "confidence intervals") values within each time bin.

```R
#Pivot the table to stack the rate estimates into a single column
coalescent <- pivot_longer(coalescent, c(growthRate1, growthRate2, growthRate3,
                                   growthRate4),
                         names_to = "time_bin", names_prefix = "growthRate")

#Summarise the rates in the log
coalescent_summary <- coalescent %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))
```

We will also add a column which converts the numbers of our time bins into their interval names: as coalescent models run from the phylogeny tips backwards in time, our "first" time bin is the youngest and our "fourth" time bin is the oldest. We will also tell R that these intervals have a set order, and specify it.

```R
#Add the interval names
coalescent_summary$interval <- c("Late Cretaceous", "Early Cretaceous",
                                 "Jurassic", "Triassic")
				 

#Ensure that the time intervals plot in the correct order
coalescent_summary$interval <- factor(coalescent_summary$interval,
                                      levels = c("Triassic", "Jurassic",
                                                 "Early Cretaceous",
                                                 "Late Cretaceous"))
```

We can plot our skylines as error bars, with a discrete bar showing the range of estimated diversification rates in each time interval.

```R
#Plot diversification skyline as error bars
ggplot(data = coalescent_summary, aes(x = interval, y = median, ymin = lowCI,
                             ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Diversification rate") +
  theme_classic(base_size = 17)
 ```
 
Alternatively, we can plot the skylines as continuous by extending our estimates across the temporal duration of each time interval, in a **piecewise constant** skyline.

```R
#Plot diversification skyline as a ribbon plot
ages <- seq.int(252, 66)
interval <- c(rep("Triassic", ((252 - 202) + 1)),
              rep("Jurassic", ((201 - 146) + 1)),
              rep("Early Cretaceous", ((145 - 101) + 1)),
              rep("Late Cretaceous", ((100 - 66) + 1)))
age_table <- as.data.frame(cbind(ages, interval))
to_plot <- left_join(coalescent_summary, age_table, by = "interval")
to_plot$ages <- as.numeric(to_plot$ages)

ggplot(to_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  xlab("Age (Ma)") + ylab("Diversification rate") +
  theme_classic(base_size = 17)
 ```
 
We can also investigate the other parameter estimated in the model, the estimated **effective population size** at the start of the coalescent process (which is the youngest end of the time interval). For us, this corresponds to an estimate of the **total species diversity** of non-avian dinosaurs at the Cretaceous-Paleogene boundary. Again, we can estimate the median and 95% HPD values for our diversity estimates.

```R
#Extract estimated diversity
pop_data <- pull(coalescent, "ePopSize")

#Summarise the diversity estimates in the log
pop_data <- as.data.frame(rbind(c(median(pop_data), quantile(pop_data, 0.025),
                                  quantile(pop_data, 0.975))))
colnames(pop_data) <- c("median", "lowCI", "highCI")
print(pop_data)
```

###The Fossilised-Birth-Death model results
We can now examine the results of our **fossilised-birth-death** model. Again, we need to read in the relevant log file and trim off the first 10% as burn-in.

```R
# Navigate to Session > Set Working Directory > Choose Directory (on RStudio)
# or change file name to the full path to the log file
#(Use "dinosaur_BDSKY_final.log" if you used our pre-cooked XML)
fbd_file <- "dinosaur_BDSKY.log"

#Read in coalescent log and trim 10% burn-in
fbd <- read.table(fbd_file, header = T) %>% slice_tail(prop = 0.9)
```

We parameterised this model using **birth**, **death** and **fossil sampling** rates, so these are the parameters included in our log. However, our exponential coalescent model estimated **diversification** rates. If we want to make the fairest comparison between the two models, we should do so using the same parameters. Fortunately, it is very straightforward to convert birth and death rates into diversification rates: simply, {% eqinline diversification = births - deaths %}. We can also calculate **turnover**, which describes the average duration of lineages (species) in our clade. Turnover is the ratio between births and deaths, so we will calculate it using {% eqinline turnover = births / deaths %}.

```R
#Calculate diversification and turnover
birth_rates <- select(fbd, starts_with("birthRate"))
death_rates <- select(fbd, starts_with("deathRate"))

div_rates <- birth_rates - death_rates
colnames(div_rates) <- paste0("divRate.",
                              seq(1:ncol(div_rates)))

TO_rates <- birth_rates / death_rates
colnames(TO_rates) <- paste0("TORate.",
                             seq(1:ncol(TO_rates)))
```

We can then calculate the median and 95% HPD values for our diversification and turnover estimates, just as we did before. This time, note that because the birth-death model runs from the origin forwards in time, our "first" time bin is the oldest and our "fourth" time bin is the youngest.

```R
#Pivot the table to stack the rate estimates into a single column
div_data <- pivot_longer(div_rates,
                         c(divRate.1, divRate.2, divRate.3, divRate.4),
                         names_to = "time_bin", names_prefix = "divRate.")

#Summarise the diversification estimates
div_data <- div_data %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))

#Add interval names
div_data$interval <- c("Triassic", "Jurassic", "Early Cretaceous",
                       "Late Cretaceous")

#Pivot the table to stack the rate estimates into a single column
turn_data <- pivot_longer(TO_rates, c(TORate.1, TORate.2, TORate.3, TORate.4),
                          names_to = "time_bin", names_prefix = "TORate.")

#Summarise the turnover estimates
turn_data <- turn_data %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))

#Add interval names
turn_data$interval <- c("Triassic", "Jurassic", "Early Cretaceous",
                        "Late Cretaceous")
```

We can also use this same process to examine our fossil sampling rates.

```R
#Pivot the table to stack the rate estimates into a single column
samp_data <- select(fbd, starts_with("sampling"))
samp_data <- pivot_longer(samp_data, c(samplingBDS.1, samplingBDS.2,
                                       samplingBDS.3, samplingBDS.4),
                          names_to = "time_bin", names_prefix = "samplingBDS.")

#Summarise log
samp_data <- samp_data %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))

#Add interval names
samp_data$interval <- c("Triassic", "Jurassic", "Early Cretaceous",
                        "Late Cretaceous")
```

Once again, we will make sure that we specify the order of our time intervals.

```R
#Ensure that the time intervals plot in the correct order
div_data$interval <- factor(div_data$interval,
                            levels = c("Triassic", "Jurassic",
                                       "Early Cretaceous", "Late Cretaceous"))

turn_data$interval <- factor(turn_data$interval,
                             levels = c("Triassic", "Jurassic",
                                        "Early Cretaceous", "Late Cretaceous"))

samp_data$interval <- factor(samp_data$interval,
                             levels = c("Triassic", "Jurassic",
                                        "Early Cretaceous", "Late Cretaceous"))
```

Again we have provided code so that the skylines can be plotted as discrete error bars...

```R
#Plot skylines as error bars
ggplot(data = div_data, aes(x = interval, y = median, ymin = lowCI,
                                      ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Diversification rate") +
  theme_classic(base_size = 17)

ggplot(data = turn_data, aes(x = interval, y = median, ymin = lowCI,
                            ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 1), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Turnover rate") +
  theme_classic(base_size = 17)

ggplot(data = samp_data, aes(x = interval, y = median, ymin = lowCI,
                            ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Sampling rate") +
  theme_classic(base_size = 17)
```

...but also as **piecewise constant** skylines using ribbon plots.

```R
#Plot skylines as a ribbon plot
ages <- seq.int(252, 66)
interval <- c(rep("Triassic", ((252 - 202) + 1)),
              rep("Jurassic", ((201 - 146) + 1)),
              rep("Early Cretaceous", ((145 - 101) + 1)),
              rep("Late Cretaceous", ((100 - 66) + 1)))
age_table <- as.data.frame(cbind(ages, interval))

div_plot <- left_join(div_data, age_table, by = "interval")
div_plot$ages <- as.numeric(div_plot$ages)
ggplot(div_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  xlab("Age (Ma)") + ylab("Diversification rate") +
  theme_classic(base_size = 17)

turn_plot <- left_join(turn_data, age_table, by = "interval")
turn_plot$ages <- as.numeric(turn_plot$ages)
ggplot(turn_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 1), colour = "black") +
  xlab("Age (Ma)") + ylab("Turnover rate") +
  theme_classic(base_size = 17)

samp_plot <- left_join(samp_data, age_table, by = "interval")
samp_plot$ages <- as.numeric(samp_plot$ages)
ggplot(samp_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  xlab("Age (Ma)") + ylab("Sampling rate") +
  theme_classic(base_size = 17)
```

We can also examine the last parameter in our fossilised-birth-death model, which is the **origin** of our clade, dinosaurs.

```R
#Extract origin data
origin_data <- pull(fbd, "origin")

#Summarise the origin estimates in the log
origin_data <- as.data.frame(rbind(c((median(origin_data) + 66),
                                     (quantile(origin_data, 0.025) + 66),
                                     (quantile(origin_data, 0.975) + 66))))
colnames(origin_data) <- c("median", "lowCI", "highCI")
print(origin_data)
```

###Model comparison
If we compare the skylines generated using our exponential coalescent and fossilised-birth-death models, we can see that they do not reconstruct the same diversification trajectory for dinosaurs across their evolutionary history. The exponential coalescent model indicates that the dinosaurs were experiencing negative diversification (net loss of species) during the Late Cretaceous, even prior to their total extinction at the Cretaceous-Paleogene boundary. In contrast, the fossilised-birth-death model suggests that the dinosaurs may have been in decline in the Early Cretaceous, but had returned to net diversification by the Late Cretaceous. One possible reason for this is that the fossilised-birth-death model includes fossil sampling rate as a parameter, whereas the exponential coalescent model simply assumes that there is no relationship between the sampling process and species richness, which may not be the case. Understanding the assumptions of the models we apply, and how well these assumptions fit our data, is therefore essential in Bayesian phylodynamics.

----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Tutorial-Template/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Tutorial-Template/master-refs.bib %}

