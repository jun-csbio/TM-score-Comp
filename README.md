![atlas](figures/logo.png)
# TM-score-Comp
TM-score-Comp: a quick and accurate method for assessing qualities of complex structure predictions of proteins, nucleic acids, and small molecule ligands

Rapid advances in AI-driven molecular structure prediction have outpaced the development of robust evaluation tools, making accurate and scalable structure comparison a critical bottleneck. Structural comparison between predicted and native complexes remains the essential gold standard, yet it becomes increasingly challenging for heterogeneous assemblies involving proteins, nucleic acids, and small-molecule ligands. Existing methods either lack a unified, size-independent scoring framework or become impractical for complexes with homologous components and symmetric ligands, where correspondence ambiguity leads to prohibitive computational cost. We present TM-score-Comp, an open-source framework for efficient and accurate evaluation of multi-component complexes. TM-score-Comp generalizes TM-score to heterogeneous assemblies and introduces two complementary metrics—rTM-score to assess molecular orientation and iTM-score to quantify interface accuracy—together enabling multi-level structural assessment. By combining orientation-aware molecule mapping with topology-guided atom mapping, TM-score-Comp achieves improved accuracy and >5× speedup over state-of-the-art methods, providing a scalable and standardized evaluation tool for complex structure prediction. The webserver of TM-score-Comp is available at https://zhanggroup.org/TM-score-Comp/.

## Installation:
* Download this repository at https://github.com/jun-csbio/TM-score-Comp.git. Then, uncompress it and run the following command lines on Linux System.

~~~
  $ cd TM-score-Comp-main
  $ chmod 777 ./make.sh
  $ ./make.sh
~~~

## Run TM-score-Comp
~~~
  $ ./exe/TMscoreC -h
  $ ./exe/TMscoreC ./example/pred.pdb ./example/pdb.pdb
~~~

## Contributing
If you encounter problems using TM-score-Comp, feel free to contact Dr. Jun Hu (hj@ism.cams.cn)! We also welcome pull requests from the community.

## References
Jun Hu, Weikang Gong, Biao Zhang and Yang Zhang. TM-score-Comp: a quick and accurate method for assessing qualities of complex structure predictions of proteins, nucleic acids, and small molecule ligands. XXXX, XX(XX): XXXX-XXXX.
