# molecular-modeling
##### Course Conclusion Paper for the Graduation Course in Computer Science from Universidade Federal de Uberlândia

----
## abstract
Programs that are considered the state of art can predict three-dimensional structures of up to 56 amino acids using force field concepts such as electrostatic attraction, van der Waals, torsion angle constraints of carbon alpha and beta chains of amino acids and bond length. 
We studied an alternative way to predict the approximate three-dimensional structure of peptides and developed a program that makes this prediction with an amino acid sequence using a Genetic Algorithm (GA). We compared the results obtained through the program with known structures, extracted from the PubChem database, comparing the total electronic energy and monitoring the execution time as we increase the size of the molecule. We get very close results in terms of energy, but the runtime increases significantly from the input of two to three amino acids. This is already expected because the problem of Molecular Modeling belongs to the NP-Complete class and the time to solve it increases exponentially as we increase the size of the molecule.

----
## usage
1. Clone the repository.
2. Install the third-party libraries.
3. Run 

----
## third-party libraries
* [Psi4NumPy](https://github.com/psi4/psi4numpy)
* [SciPy](https://github.com/scipy/scipy)
