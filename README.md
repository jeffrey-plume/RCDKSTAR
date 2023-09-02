RCDKSTAR
================

## Introduction

**RCDKSTAR** stands for R-Chemistry Development Kit Structure-Activity
Relationships. RCDKSTAR is a lightweight BioAssay aggregator and
structural-activity relationship modeling suite.

## Installation

``` r
shiny::runGitHub( "RCDKSTAR", "jeffrey-plume")
```

## Modeling with RCDKSTAR

Select a protein coding gene from the drop-down list to populate the
screen with summary data from PubChem BioAssays targeting the protein
product. Filter the data as desired. Error and warning messages are
displayed if something is off. BioAssays with the same target may have
different endpoints, which can affect the activity outcome of the
experiments. Choose a method for fingerprinting the molecules and any
additional calculated properties. Hold CTRL to select multiple
additional properties. Press ‘Load Properties’ to retrieve the desired
values. Preprocessing of the data can be done by checking any of the
checkbox options. Once satisfied, click “Train Model”. Depending on the
complexity of the model and computing power, this may take a while. By
Default, 75% of the curated dataset is used for training and the
remainder for validation. Once complete, tables displaying summary
statistics, rank-order motifs by importance, and probability outcomes of
the remaining 25% of the data not used to train the model. Click “Test
Unknowns” to apply the model to other compounds. If no SDF compound
structural file is imported, compounds of the Drug Repurposing Hub by
the Broad Institute are analyzed for activity.

<div align="center">

<img src="Peek 2022-06-29 19-28.gif" width="67%" />

</div>
