![atlas](figures/logo.png)
# TM-score-Comp
TM-score-Comp: a quick and accurate algorithm for measuring quality of complex structure predictions of proteins, nucleic acids, and small molecule ligands

Structural comparison between native and predicted structures is the gold-standard way for accessing the quality of 3D structure predictions. As structure prediction advances to include complex structures with different molecule types, developing structure comparison methods to accurately assess predictions for diverse molecules (proteins, nucleic acids, small molecules) in various forms (monomers, oligomers, and polymers) is critical for advancing structural bioinformatics. However, the existing structure comparison methods cannot directly evaluate predicted complexes consisting of proteins, nucleic acids, and small molecule ligands, nor can they provide optimal/suboptimal molecule mappings for homologous polymers within a limited timeframe. We have developed a new and open-source structure comparison method, TM-score-Comp, to quickly and accurately measure the quality of structure predictions of complexes composed of one or more molecule types, such as proteins, DNAs, RNAs, and small molecule ligands. Large-scale benchmarks demonstrated consistent advantages of TM-score-Comp over state-of-the-art methods in the assessment of 3D structure predictions of different molecule types and runs more than 5 times faster on average, especially for those with homologous molecules. Detailed analyses demonstrated that the main advantage of TM-score-Comp lies in the molecule mapping algorithm can always give an optimal/suboptimal molecule mapping results in fast, resulting the improvement of the accuracy and speed of the complex structure comparison process. The on-line web server and source code of TM-score-Comp are made freely available at <b>https://zhanglab.comp.nus.edu.sg/TM-score-Comp/</b> for academical use.

## Installation:
* Download this repository at https://github.com/jun-csbio/TM-score-Comp.git. Then, uncompress it and run the following command lines on Linux System.

~~~
  $ cd TM-score-Comp-main
  $ chmod 777 ./make.sh
  $ ./make.sh
~~~

* If you want run TM-score-Comp in Windows System. You can directly run the files of "./exe/TMscoreC.exe" or "./exe/TMscoreCmt.exe".

## Run Single-Thread TM-score-Comp
~~~
  In Linux System
  $ ./exe/TMscoreC -h
  $ ./exe/TMscoreC ./example/pdb.pdb ./example/pred.pdb

  or in Windows System
  $ ./exe/TMscoreC.exe -h
  $ ./exe/TMscoreC.exe ./example/pdb.pdb ./example/pred.pdb
~~~

## Run Multi-Thread TM-score-Comp
~~~
  In Linux System
  $ ./exe/TMscoreCmt -h
  $ ./exe/TMscoreCmt ./example/pdb.pdb ./example/pred.pdb

  or in Windows System
  $ ./exe/TMscoreCmt.exe -h
  $ ./exe/TMscoreCmt.exe ./example/pdb.pdb ./example/pred.pdb
~~~

## Example of Result
![atlas](figures/output_example.png)

## Contributing
If you encounter problems using TM-score-Comp, feel free to create an issue! We also welcome pull requests from the community.

## References
Jun Hu, Weikang Gong, Biao Zhang and Yang Zhang. TM-score-Comp: a quick and accurate tool for assessing quality of complex structure predictions of proteins, nucleic acids, and small molecule ligands. XXXX, XX(XX): XXXX-XXXX.
