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

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### Tracer

[Tracer] (http://beast.community/tracer) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v1.7.x.

### R / RStudio

We will be using R to analyze the output of the Birth-Death Skyline plot. RStudio provides a user-friendly graphical user interface to R that makes it easier to edit and run scripts. (It is not necessary to use RStudio for this tutorial).

----

# Practical: Skyline analyses for macroevolution

In this tutorial we will estimate diversification rates for dinosaurs using a previously published supertree.

The aim of this tutorial is to:
- Learn how to set up a skyline analysis using a previously made phylogeny;
- Develop skills in simple xml hacking;
- Highlight the differences between exponential coalescent and fossilised-birth-death skylines.

## The data
We will be inferring our skyline parameters using a ready-made phylogeny containing 420 dinosaur species, published by {% cite Lloyd2008 --file Tutorial-Template/master-refs.bib %}. This phylogeny is a 'supertree', created using an informal method to collate several smaller dinosaur phylogenies into a larger one. Supertrees are typically cladograms, and must be timescaled using fossil data to estimate their branch lengths. The branch lengths we use here were inferred by {% cite Sakamoto2016 --file Tutorial-Template/master-refs.bib %}: they used fossil occurrences from the [Paleobiology Database] (http://paleobiodb.org) to infer the midpoint of the temporal range of each tip, and timescaled the phylogeny using the "equal" method, which distributes time evenly between the available branches.

## Creating an xml template with BEAUti
Many of the features we will need in our xml files are currently not implemented in BEAUti. However, we will start our analyses by creating an xml in BEAUti which will then serve as a template for us to alter by hand ("hack") later.

### Install BEAST2 packages
The coalescent-based skyline model is included in the core of BEAST2, but we need to install the **BDSKY** package, which contains the birth-death skyline model. We will also need the **feast** package, which will allow us to integrate some more complex features into our xml. Installation of packages is done using the package manager, which is integrated into BEAUti.

>Open the **BEAST2 Package Manager** by navigating to **File > Manage Packages**.
>
>Install the **BDSKY** and **feast** packages by selecting them and clicking the **Install/Upgrade** button one at a time.

After the installation of a package, the program is on your computer, but BEAUti is unable to load the template files for the newly installed packages unless it is restarted.

>Close the **BEAST2 Package Manager** and **_restart_** BEAUti to fully load the **BDSKY** and **feast** packages.

### Setting up the Fossilised-Birth-Death skyline analysis
The first step in setting up an analysis in BEAUti is to upload an alignment. In our case, we do not have an alightment as we are using a ready-made phylogeny instead, but we will still need to upload an alignment in order to initialise our xml. We have provided a dummy `.nexus` file for this purpose, which we will replace later with code that reads in our phylogeny.

>In the **Partitions** panel, import the nexus file with the empty alignment by navigating to **File > Import Alignment** in the menu and then finding the `empty.nexus` file on your computer, *or* drag and drop the file into the **BEAUti** window.

The empty alignment should appear in the **Partitions** tab, named `empty` and containing a single nucleotide site for a single taxon.

Most of the default tabs in BEAUti relate to inferring the phylogeny, which we will not be doing, so we can skip straight to the **Priors** tab. It's worth remembering that a lot of the default parameters in BEAUti are set up for analyses very different to ours, so it's worth checking these thoroughly.

>Select the **Priors** tab and choose **Birth Death Skyline BDSParam** as the tree prior.

The "standard" `contemporary` birth death skyline is parameterised using a reproductive number, a "become uninfectious" rate and a sampling proportion. However, here we are using the `BDSParam` model, which refers to parameterisation using a birth rate, a death rate and an extant sampling proportion. As our phylogeny is for non-avian dinosaurs, for which we only have fossils and no extant (genetic) samples, we will later exchange the extant sampling proportion, usually denoted using {% eqinline rho %}, for a fossil sampling rate, usually denoted using {% eqinline psi %}.

Choosing sensible priors for these parameters is not straightforward, so we will select priors which are relatively unrestrictive. For the birth and death rates, we will use exponential priors with a mean of 1.0; this places more probability on small rates, but still permits rates which are towards the higher end of those estimated from living animals and plants {% cite HenaoDiaz2019 --file Tutorial-Template/master-refs.bib %}. For the sampling rate, we will also choose an exponential prior, this time with a mean of 0.2.

On the **initial =** buttons you will see two sets of square brackets. The first indicates what the starting value for that parameter will be, meaning its value in the first iteration of your chain. We need to alter our initialisation values to ensure that they sit within our prior distributions. The second contains two values which denote the range of values that our parameters are permitted to take, which are also important to consider carefully. Priors influence the probability of certain values being tested in your chain, but values which are improbable under your prior can still be selected if the signal in your data is strong enough; setting this range provides hard limits to your parameter values regardless of your priors. Our parameters are all rates, expressed per branch per million years, so the full set of values they can take ranges between 0.0 and infinity.

This is also the point where we express how many sections (here called **dimensions**) we want in our **piecewise constant** rates. Our rates will be assumed to be constant within these sections, but will be permitted to change at the break points between them. To keep our analysis simple (in the hope of a timely convergence!), we are going to give our rates **four** dimensions, corresponding to the major geological intervals spanned by the phylogeny: the Triassic, Jurassic, Early Cretaceous and Late Cretaceous.

>Change the **birth rate** prior from a uniform distribution to an exponential. Using the drop-down arrow on the left, check that the **mean** is set to 1.0. Click on the **initial =** button and change the **initialisation value** to 1.0. Check that the lower value is 0.0 and the upper value is `Infinity`. Change the **Dimension** to 4. Repeat these four steps for the **death rate** prior.
>
>Change the **rho** (sampling) prior from a uniform distribution to an exponential. Using the drop-down arrow on the left, change the **mean** value to 0.2. Click on the **initial =** button and change the **initialisation value** to 0.5. Check that the lower value is 0.0, and change the upper value to `Infinity`. Change the **Dimension** to 4.



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

