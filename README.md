# Mass Decomposition
This algorithm performs the decomposition of a positive real-numbered value over a weighted alphabet.

There are two parts to this algorithm in it's current ideation:
* The creation of the universal search space.
* The querying of said search space without duplicates using a branch and bound algorithm that reduces the size of the search space.

## Usage
The standalone algorithm implementation reguires the target value (query mass) and a single PTM as the only arguments. If there is no PTM, then that argument should be the character 'x'.
The output in the provided implementation is written to the terminal. If desired, this can be changed to write to file which would make sense with larger target values that produce many millions of sequences.
The executable is compiled using the standard gcc compiler, using the -lm flag to include math.h successfully.

&copy; Christopher Burgoyne 2021
