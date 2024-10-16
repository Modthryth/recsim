#include "altRandom.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <fstream>
//Recombination: two distinct oligomers chosen at random and recombine at a random bond
//AB + CD = AC + BD
//For example, 6 + 4 can be 9/1, 8/2, 7/3, 6/4 (yes the same), or 5/5, reverse order is redundant for size
//In essence, first oligo is random of (AB + CD) and next oligo makes up the difference
//The first chosen oligo (bigger or smaller) either decays at a random bond or recombines
//p = probability of a nucleotide adding to the chain
//No = initial number of monomers
//N = total number of polymers (including monomers)
//N(avg) = 1 / (1 - p)
//p = (No + N) / No
//pL = probability of chain length L = p ^ L
//NpL = number of polymers of size L = N*p^L

using std::cout;
using std::endl;
using std::cin;
using std::vector;
using std::sort;

//Using Arvid Gestmann's random number generator algorithm PCG for speed
class PCGRandomGenerator { //To implement any of altRandom's generators - simply replace pcg eng with xorshift eng, etc
public:
	PCGRandomGenerator() : eng(rd) {}

	//Gets a random double
	double GetRandomDbl() {
		std::uniform_real_distribution<double> dDistribution(0, 1);
		return dDistribution(eng);
	}

	//Gets a random int, max is exclusive upper bound
	int GetRandomInt(int min, int max) {
		std::uniform_int_distribution<int> iDistribution(min, max - 1);
		return iDistribution(eng);
	}

	//Gets a random int not equal to iAvoid
	int GetRandomInt(int min, int max, int iAvoid) {
		int returnVal = 1;
		do {
			std::uniform_int_distribution<int> iDistribution(min, max - 1);
			returnVal = iDistribution(eng);
		} while (returnVal == iAvoid);
		return returnVal;
	}
private:
	std::random_device rd;
	pcg eng; //choice of engine
};

double distVal = 100000000.0 / 1048575.0; //exponential starts with 100 million total oligomers. distVal is equal to about 95 and represents the amount of 20mers which is the maximum size allowed to start 

int RecombinableSize = 1; //minimum recombinable size (MRS)
int mintargetcomp = -1; //Minimum array index for the MRS. Random oligomers are chosen from all length classes with an index greater than this value (index 0 = length 1, or MRS - 2)

bool gradedRec = false; //If smaller oligos have a lower probability of recombination instead of total exclusion. Requires MRS = 1
double gradedRecProbabilities[] = { 0.2, 0.4, 0.6, 0.8 }; //For sizes 1-4

int outX = 5; //the number of replicates

int AvailableOligomerSum = 0; //The total number of oligomers available to recombine (must be equal to or greater than RecombinableSize)
double initpHyd = 0.0; //Rate of hydrolysis out of 1.0, the remaining difference results in recombination. Ranges from 0.2 to 0.8 in published simulations
int cyclemax = 500000000; //Total number of cycles to run if the pool is not expected to decay. Limited by integer max , approx 2.1 x 10^9

bool writeDist = true; //whether to output the distribution, takes time and generates a lot of data
bool writeMax = false; //whether to write the tracked maximum sizes to output files

//polymerization variables
bool polymerize = false; //whether to do single nucleotide polymerization 

double p = 0.75; //probability of a nucleotide adding to the chain; average length = 1 / (1 - p)
int Imonomers = 100000000; //the number of monomers to start (No)
int totalpolymers = (Imonomers * (1.0 - p)); //number of polymers N; p = (No - N) / No
int Nmonomers = totalpolymers * (1.0 - p); //number of expected monomers after polymerization

int polyRate = 1; //polymerization rate = number of rec events that occur per one polymerization event. Set to 1 for highest rate, larger numbers for lower rate. Cannot be 0.
bool altPolyStart = true; //use polymerization algorithm to generate distribution instead of fixed
bool fixed_polymer_start = false; //whether to use a selection of single-length polymers to start - for example only 10mers

int num100mers = 0; //how many 100mers or higher are formed by recombination

std::ofstream opfile; //Output file for writing to comma delimited
std::ofstream distfile; //file for writing distribution data
std::ofstream maxfile; //for writing maximum size;

int main(int argc, char* argv[]) {

	PCGRandomGenerator Rand;
	int distribution[1000]; //max size up to 1000 nt, but track
	int saved_dist[1000]; //carbon copy for replenishment purposes
	double avgdist[1000]; //for computing averages
	double dvals[20]; //Individual oligo probabilities for initial distribution

	//Initialization of variables
	for (int i = 0; i < 20; i++) { dvals[i] = 0.0; }
	int incr = 1;
	int distnum = (int)(distVal + 0.5);

	for (int i = 0; i < 1000; i++) {
		distribution[i] = 0;
		saved_dist[i] = 0;
		avgdist[i] = 0.0;
	}

	//Open files
	opfile.open("RecSimData.txt");
	opfile.close();

	distfile.open("RecSimDist.txt");
	distfile.close();

	if (writeMax) {
		maxfile.open("MaxSize.txt");
		maxfile.close();
	}

	for (int suploop = 0; suploop < 1; suploop++) { //super loop, coordinates with pHyd to do sequential simulations at different values
		double pHyd = initpHyd - suploop * 0.1; //coordinate with super loop
		double avgmax = 0.0; //average maximum size
		double cycleavg = 0.0; //Average number of cycles until decay
		int globalmax = 0; //The largest size obtained in all replicates

						   //Write file header
		opfile.open("RecSimData.txt", std::ios_base::app);
		opfile << "Hydrolysis rate: " << pHyd << endl;
		opfile << "Maximum allowed oligo: " << RecombinableSize << endl;
		opfile << "Maximum Size, Events" << endl;
		opfile.close();

		if (writeMax) {
			maxfile.open("MaxSize.txt", std::ios_base::app);
			maxfile << "Hydrolysis rate: " << pHyd << endl;
			maxfile << "Maximum allowed oligo: " << RecombinableSize << endl;
			maxfile << "Maximum Size, Events" << endl;
			maxfile.close();
		}

		//Outer loop for outX replicates, set to 5 in published simulations
		for (int oloop = 0; oloop < outX; oloop++) {

			//Variables 
			int RecEvents = 0;
			int maxsize = 20;
			int maxloopcount = 0;
			num100mers = 0;
			int loopsum = 0;

			//reset distribution
			int incr = 1;
			int distnum = (int)(distVal + 0.5);

			for (int i = 0; i < 1000; i++) {
				distribution[i] = 0;
			}

			//Set up the starting distribution
			if (altPolyStart) {
				//Polymerize the distribution using the polymerase algorithm
				distribution[0] = 100000000;
				do {
					int clen = 1;
					distribution[0]--; //remove initial which is now clen
					double gp = 0.0;
					do {
						if (distribution[0] > 0) {
							gp = Rand.GetRandomDbl();
							if (gp < p) { clen++; distribution[0]--; }
						}
					} while (gp < p);
					distribution[clen - 1]++;
					if (clen > maxsize) { maxsize = clen; } //track max size
				} while (distribution[0] > Nmonomers);
			}
			else {
				//Fixed length distribution (assigned)
				if (fixed_polymer_start) {
					distribution[4] = 20000000;
					saved_dist[4] = 20000000;
				}
				else {
					//Fixed exponential distribution
					for (int i = 19; i > -1; i--) {
						distribution[i] = distnum;
						saved_dist[i] = distnum;
						dvals[i] = (double)distnum / 100000000.0;

						incr *= 2;
						distnum = (int)(distVal * incr + 0.5); //Iterate up instead of down using distVal, so that we reach approximately 50,000,000 monomers
					}
				}
			}

			//Write initial distribution (zero timepoint)
			if (writeDist) {
				distfile.open("RecSimDist.txt", std::ios_base::app);
				distfile << "Next replicate: " << oloop << endl;
				for (int ibc = 0; ibc < 200; ibc++) {
					distfile << ibc << ",";
				}
				distfile << endl;
				distfile << RecEvents << ",";
				for (int ibc = 0; ibc < maxsize; ibc++) {
					int outw = distribution[ibc];
					distfile << outw << ",";
				}
				distfile << endl;
				distfile.close();
			}

			//inner loop
			//do until nothing left or 300,000,000 events
			do {
				AvailableOligomerSum = 0;
				//Choose at random
				//Get total of recombinable oligos
				for (int i = 0; i < 1000; i++) {
					if (i > mintargetcomp) { AvailableOligomerSum += distribution[i]; } // AvailableOligomerSum is the total of all recombinable oligos
				}
				loopsum = AvailableOligomerSum;
				if (AvailableOligomerSum > 1) { //have to have at least two oligos
					int inext = Rand.GetRandomInt(0, AvailableOligomerSum); //randomly choose
					int targlen = 1;
					for (int i = 0; i < 1000; i++) {
						if (i > mintargetcomp) { inext -= distribution[i]; } //Iterate over distribution to find the length of the chosen oligo
						if (inext < 0) { targlen = i + 1; break; }
					}
					//temporarily remove the first one to avoid influencing probability of selection
					distribution[targlen - 1] -= 1;

					AvailableOligomerSum = 0; //reset
											  //Now get second oligo by same manner
					for (int i = 0; i < 1000; i++) {
						if (i > mintargetcomp) { AvailableOligomerSum += distribution[i]; }
					}
					int inext2 = Rand.GetRandomInt(0, AvailableOligomerSum);
					int targlen2 = 1;
					for (int i = 0; i < 1000; i++) {
						if (i > mintargetcomp) { inext2 -= distribution[i]; } //get size of second oligomer
						if (inext2 < 0) { targlen2 = i + 1; break; }
					}
					//Remove the second oligo
					distribution[targlen2 - 1] -= 1;

					//Recombine or break oligos, then add back the one removed
					double dprob = Rand.GetRandomDbl();
					if (dprob < pHyd) { //break probability must be between 0 and 1
										//Break only
										//Breakage requires minimum size of two, or one bond, a monomer can't break
						if (targlen > 1) {
							int ires = Rand.GetRandomInt(1, targlen);
							distribution[ires - 1]++;
							distribution[targlen - ires - 1]++;
							distribution[targlen2 - 1] += 1; //If first gets hydrolyzed, add the second back to the distribution since nothing happened
						}
						else {
							distribution[targlen - 1] += 1; //Add back since a monomer can't be hydrolyzed
							distribution[targlen2 - 1] += 1; //Add second oligo back 
						}
					}
					else {
						//Break and join (two monomers has no effect since ires will be 1)
						int ires = Rand.GetRandomInt(1, (targlen + targlen2)); // Get random location of their combined lengths

						//Check recombinable size
						bool doRec = (targlen >= RecombinableSize) && (targlen2 >= RecombinableSize);

						if (gradedRec) {
							//Graded recombination
							//Reduced probability if either of the two oligos is below 5 nt, with the probability defined by the smaller one
							double gradprob = 1.0;
							int minsize = targlen > targlen2 ? targlen2 : targlen;
							if (minsize < 5) { gradprob = gradedRecProbabilities[minsize - 1]; }
							double grade = Rand.GetRandomDbl();
							doRec = doRec && (grade < gradprob);
						}

						//Restrict to 1000 max
						if ((ires < 1000) && (doRec)) {
							if (ires > maxsize) { maxsize = ires; } //Track largest oligomer
							if (ires >= 100) { num100mers++; }
							distribution[ires - 1]++;
							distribution[targlen + targlen2 - ires - 1]++;
						}
						else {
							//won't allow if too big for array purposes or below MRS, so add back originals
							distribution[targlen - 1] += 1;
							distribution[targlen2 - 1] += 1;
						}
					}
				}
				else {
					//if no more oligos left, then done
					cycleavg += (double)RecEvents; //how many events it took to finish the pool
					break;
				}
				RecEvents++;

				//Ppolymerization
				if (polymerize) {
					if (RecEvents % polyRate == 0) {
						int clen = 1;
						double gperc = distribution[0] / (double)Imonomers;
						double cval = Nmonomers / (double)Imonomers;
						if ((distribution[0] > 1) && (gperc >= cval)) { //has to have the minimum concentration
							distribution[0]--; //remove initial which is now clen
							double gp = 0.0;
							do {
								if (distribution[0] > 0) { //bounds checking
									gp = Rand.GetRandomDbl();
									if (gp < p) { clen++; distribution[0]--; } //add to the chain if it meets the probability
								}
							} while (gp < p);
							if ((clen) > maxsize) { maxsize = clen; } //track largest
							distribution[clen - 1]++;
						}
					}
				}

				if (writeDist) {
					if ((RecEvents % 1000000) == 0) {
						//Write distribution every 100k events (change to larger number if minimum size is smaller to save time and memory)
						distfile.open("RecSimDist.txt", std::ios_base::app);
						distfile << RecEvents << ",";
						for (int ibc = 0; ibc < maxsize; ibc++) {
							int outw = distribution[ibc];
							distfile << outw << ",";
						}
						distfile << endl;
						distfile.close();
					}
				}

				if (writeMax) {
					//Output the maximum size every 1 million iterations
					if ((RecEvents % 1000000) == 0) {
						maxfile.open("MaxSize.txt", std::ios_base::app);
						//get max index
						int maxs = 0;
						for (int im = 999; im > -1; im--) {
							if (distribution[im] > 0) {
								maxs = im; break;
							}
						}
						maxfile << maxs << endl;
						maxfile.close();
					}
				}

				if ((RecEvents % 10000000) == 0) { cout << "Reaction events: " << RecEvents << endl; } //keeps track while running simulations for long time periods

				if ((RecEvents % 100000000) == 0) {
					cout << "Current distribution: " << endl;
					for (int i = 0; i < 200; i++) {
						if (distribution[i] > 0) {
							cout << i + 1 << "mers: " << distribution[i] << endl;
						}
					}
				}
				maxloopcount++;
			} while ((loopsum > 1) && (maxloopcount < cyclemax)); //Loop until nothing left to recombine or hit maximum cycles

			cout << endl;
			cout << "Replicate: " << (oloop + 1) << endl;
			cout << "The maximum attained size is: " << maxsize << endl;
			cout << "The number of events until complete decay is: " << RecEvents << endl;
			for (int i = 0; i < 200; i++) {
				if (distribution[i] > 0) {
					cout << i + 1 << "mers: " << distribution[i] << endl;
				}
			}
			opfile.open("RecSimData.txt", std::ios_base::app);
			opfile << maxsize << "," << RecEvents << endl;
			opfile.close();
			avgmax += (double)maxsize;

			//Write final distribution
			if (writeDist) {
				distfile.open("RecSimDist.txt", std::ios_base::app);
				distfile << RecEvents << ",";
				for (int ibc = 0; ibc < maxsize; ibc++) {
					int outw = distribution[ibc];
					//if (distribution[ibc] < 1) { outw = 1; }
					distfile << outw << ",";
				}
				distfile << endl;
				//distfile << "Number of 100mers: " << num100mers << endl;
				distfile.close();
			}

			//Add final distribution to the average file
			for (int iad = 0; iad < 1000; iad++) {
				avgdist[iad] = avgdist[iad] + distribution[iad] * 1.0;
			}

			//Get global max size
			if (maxsize > globalmax) { globalmax = maxsize; }
		}

		//Compute and write out average
		for (int iad = 0; iad < 1000; iad++) {
			avgdist[iad] = avgdist[iad] / (double)outX;
		}

		if (writeDist) {
			distfile.open("RecSimDist.txt", std::ios_base::app);
			distfile << "Average Distribution of " << outX << " replicates: " << endl;
			for (int ibc = 0; ibc < globalmax; ibc++) {
				double dw = avgdist[ibc];
				if ((dw < 1.0) && (dw > 0.0)) { dw = 1.0; }
				distfile << std::fixed << dw << ",";
			}
			distfile << endl;
			distfile.close();
		}

		//rewrite max size lines here
		if (writeMax) {
			maxfile.open("MaxSize.txt", std::ios_base::app);
			maxfile << "End of iteration" << endl;
			maxfile.close();
		}

		//Re-zero the average array
		for (int iad = 0; iad < 1000; iad++) {
			avgdist[iad] = 0.0;
		}
		avgmax = avgmax / (double)outX;
		cycleavg = cycleavg / (double)(outX);
		cout << "Final statistics: " << endl;
		cout << "Replicates: " << outX << endl;
		cout << "Hydrolysis rate: " << pHyd << endl;
		cout << "Minimum reactive oligomer size: " << (mintargetcomp + 2) << endl;
		cout << "Average maximum oligomer: " << avgmax << endl;
		cout << "Average events until decay: " << cycleavg << endl;
	}
	cin.ignore();
}