This is a collection of code verifying the results of the manuscript ``"Modular Curves of Prime-Power Level with Infinitely Many Quadratic Points''``, by Michael Cerchia and Rakvi.

Files needed to run the code : Please download the ``LMFDB Data folder``. We use the FindModelofXG function available at https://github.com/davidzywina/OpenImage. We also need the congruence subgroup data taken from Cummins-Pauli database of congruence subgroups available at https://math-sites.uncg.edu/sites/pauli/congruence/. We provide the instructions to run these files as needed in individual files.

The folder ``LMFDB Data`` contains data downloaded from LMFDB database of modular curves.

The folder ``Magma Code`` contains a file for each positive rank bielliptic modular curve for which Magma is used in showing it is in fact positive rank bielliptic. The file names correspond to the LMFDB labels. 

The folder ``Not positive rank`` contains a file for each candidate curve that is not hyperelliptic and fails to be positive rank bielliptic. In the argument summaries at the top of each file, we write "not bielliptic" as shorthand for "not positive rank bielliptic". 

The file ``Bielliptic prime power level upper bound on GL2 level`` contains verifications of section 4.3.2 of paper.

The files ``Genus 0`` and ``Genus 1`` contain verifications for sections 4.1 and 4.2 respectively.

The files ``Hyperelliptic prime power level upper bound on GL2 level``, ``Hyperellipticcandidates`` and ``Remaining cases-hyperelliptic`` contain verifications for section 4.3.1 and 5.

The file ``Section 7.1`` contains verifications for section 7.1.

The file ``Table 6`` contains the code to produce table 6.

The file ``biellipticlabels twist`` contains the code to produce candidates to check for positive rank bielliptic case and also produces rank 0 curves.



