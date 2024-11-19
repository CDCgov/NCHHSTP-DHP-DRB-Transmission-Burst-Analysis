# HIV Transmission Burst Analysis

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

# Overview

This analytical approach was developed to identify bursts of rapid HIV transmission across a clock-tree and characterize the contribution of bursts of rapid transmission to future transmission. The application of this approach is described in detail in a manuscript currently under review for publication.

## Outline

The script infers a maximum likelihood phylogeny using FastTree2 and optimizes time-scaled phylogenies and infers transmission dates corresponding to dated internal nodes within each phylogeny using TreeTime. HIV transmission bursts are identified as three or more internal nodes descended from a single lineage within a defined detection period and are detected using tree-traversal implemented in ETE3. Lineages during the burst detection period are classified as burst-involved or non-burst-involved, where a lineage is considered involved in a transmission burst if it is the product of an inferred transmission event within a burst. Inferred transmission events, or internal nodes, dated during a defined follow-up period are identified as descended from lineages involved or not involved in transmission bursts during the burst detection period. 

To estimate the relative contribution of transmission bursts to future transmission, divide the number of internal nodes during the follow-up period descended from bursts (A1) by the number of lineages active during the burst detection period involved in bursts (B1). Then, compare this ratio with the number of non-burst-descended inferred transmission events during the follow-up period (A2) divided by the number of non-burst-involved lineages during the burst detection period (B2). The relative contribution ratio is estimated as: 
(Contribution of lineages involved in bursts to future transmission)/(Contribution of lineages not involved in bursts to future transmission)=(A_1/B_1)/(A_2/B_2)

To describe populations affected by transmission bursts throughout the tree, rather than bursts identified within a single, defined burst detection period, bursts may be identified using a sliding window of a defined width. Phylogenetic tips are then classified as either (i) a member or descendant of a transmission burst or (ii) not a member or descendant of a transmission burst. 

## Required Dependencies

* ete3: http://etetoolkit.org/download/
* FastTree: http://www.microbesonline.org/fasttree/#Install
* TreeTime: https://treetime.readthedocs.io/en/latest/installation.html

## Defined Parameters

* epochStart: Start year for fixed burst detection period (integer)
* epochEnd: End year for fixed burst detection period (integer)
* events: Number of inferred internal nodes descended from a single lineage during the burst detection required to identify a burst
* upperDxLimit: End year for fixed follow-up period (integer)

## Input File Requirements

* Tab-delimited file
* Unique sequence- or person-level identifier
* Cleaned subtype B HIV sequence including reverse transcriptase region
* Genotype sample date (YYYYMMDD)
* Diagnosis date (YYYYMMDD)

## Analysis Script

https://github.com/CDCgov/NCHHSTP-DHP-DRB-Transmission-Burst-Analysis/blob/main/src/runBurstAnalysis.py

# Standard Notices
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Related documents
* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
