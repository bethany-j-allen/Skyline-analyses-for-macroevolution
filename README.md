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

## Creating an XML template with BEAUti
Many of the features we will need in our XML files are currently not implemented in BEAUti. However, we will start our analyses by creating XML files in BEAUti which will then serve as a template for us to alter by hand ("hack") later.

### Install BEAST2 packages
The coalescent-based skyline model is included in the core of BEAST2, but we need to install the **BDSKY** package, which contains the birth-death skyline model. We will also need the **feast** package, which will allow us to integrate some more complex features into our analyses. Installation of packages is done using the package manager, which is integrated into BEAUti.

>Open the **BEAST2 Package Manager** by navigating to **File > Manage Packages**.
>
>Install the **BDSKY** and **feast** packages by selecting them and clicking the **Install/Upgrade** button one at a time.

After the installation of a package, the program is on your computer, but BEAUti is unable to load the template files for the newly installed packages unless it is restarted.

>Close the **BEAST2 Package Manager** and **_restart_** BEAUti to fully load the **BDSKY** and **feast** packages.

### Setting up the Fossilised-Birth-Death skyline analysis
The first step in setting up an analysis in BEAUti is to upload an alignment. In our case, we do not have an alignment as we are using a ready-made phylogeny instead, but we will still need to upload an alignment in order to initialise our XML. We have provided a dummy `.nexus` file for this purpose, which we will replace later with code that reads in our phylogeny.

>In the **Partitions** panel, import the nexus file with the empty alignment by navigating to **File > Import Alignment** in the menu and then finding the `empty.nexus` file on your computer, *or* drag and drop the file into the **BEAUti** window.

The empty alignment should appear in the **Partitions** tab, named `empty` and containing a single nucleotide site for a single taxon.

Most of the default tabs in BEAUti relate to inferring the phylogeny, which we will not be doing, so we can skip straight to the **Priors** tab. It's worth remembering that a lot of the default parameters in BEAUti are set up for analyses very different to ours, so it's worth checking these thoroughly.

>Select the **Priors** tab and choose **Birth Death Skyline BDSParam** as the tree prior.

The "standard" `contemporary` birth death skyline is parameterised using a reproductive number, a "become uninfectious" rate and a sampling proportion. However, here we are using the `BDSParam` model, which refers to parameterisation using a birth rate, a death rate and an extant sampling proportion. As our phylogeny is for non-avian dinosaurs, for which we only have fossils and no extant (genetic) samples, we will later exchange the extant sampling proportion, usually denoted using {% eqinline rho %}, for a fossil sampling rate, usually denoted using {% eqinline psi %}.

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

The first line sets out information about the format of the xml and information about its contents, which we can ignore. The next section is labelled `data`, and contains our dummy alignment. We don't need this section, so it can be commented out or simply deleted.

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

Following this we see a description of our three inferred rates: birth, death and sampling. These lines include a lot of information, much of which we specified earlier via the **Priors** tab in BEAUti. As previously mentioned, our sampling parameter is currently set up as the extant sampling proportion, `rho`, which we want to exchange for a fossil sampling rate, `sampling`. For the sake of transparency, we will change the `id` (name) of our sampling parameter here from `rhoBDS` to `samplingBDS`. Here, this is simply a label and does not change what the parameter does, but it is how the parameter is referenced in other parts of the XML, where its name will also need to be changed.

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

