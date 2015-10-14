/*
 * Copyright 2013-{year} Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.jsres;

import java.util.Arrays;
import java.util.Random;

public final class SRES {
    public static final Random random = new Random(51102);

    private static final int boundingTries = 10;
    private final Objective objective;
    private final double tau, tauDash, rankingPenalizationFactor;
    private final int lambda, mu, numberOfFeatures, numberOfSweeps;

    public SRES(Objective objective) {
        this(objective, 200, 30, 1.0, 200, 0.45);
    }

    public SRES(Objective objective,
                int lambda, int mu,
                double expectedConvergenceRate,
                int numberOfSweeps, double rankingPenalizationFactor) {
        this.objective = objective;
        this.numberOfFeatures = objective.getNumberOfFeatures();
        this.lambda = lambda;
        this.mu = mu;
        this.tau = expectedConvergenceRate / Math.sqrt(2 * Math.sqrt(numberOfFeatures));
        this.tauDash = expectedConvergenceRate / Math.sqrt(2 * numberOfFeatures);
        this.numberOfSweeps = numberOfSweeps;
        this.rankingPenalizationFactor = rankingPenalizationFactor;
    }

    public SolutionSet run() {
        return run(1750);
    }

    public SolutionSet run(int numberOfGenerations) {
        SolutionSet solutionSet = new SolutionSet();

        for (int i = 0; i < numberOfGenerations; i++) {
            solutionSet.evaluate();
            solutionSet.sort();
            solutionSet = solutionSet.evolve();
        }

        solutionSet.evaluate();
        solutionSet.sort();

        return solutionSet;
    }

    public class SolutionSet {
        private final Solution[] solutions;

        SolutionSet() {
            this.solutions = new Solution[lambda];
            for (int i = 0; i < lambda; i++) {
                solutions[i] = new Solution(objective.generateFeatureVector(), objective.getMutationRates());
            }
        }

        SolutionSet(Solution[] solutions) {
            this.solutions = solutions;
        }

        void sort() {
            for (int i = 0; i < numberOfSweeps; i++) {
                boolean swapped = false;
                for (int j = 0; j < lambda - 1; j++) {
                    double p = random.nextDouble(),
                            p1 = solutions[j].getPenalty(),
                            p2 = solutions[j + 1].getPenalty();
                    if (p < rankingPenalizationFactor || (p1 == 0 && p2 == 0)) {
                        // no solution break constraints or penalization is not applied
                        if (solutions[j].getFitness() > solutions[j + 1].getFitness()) {
                            swap(j, j + 1);
                            swapped = true;
                        }
                    } else if (p1 > p2) {
                        // constraint penalization
                        swap(j, j + 1);
                        swapped = true;
                    }
                }
                if (!swapped) {
                    // sorted
                    break;
                }
            }
        }

        SolutionSet evolve() {
            Solution[] newSolutions = new Solution[lambda];

            int j = 0;
            for (int i = 0; i < lambda; i++) {
                // mutation rates are sampled from top solutions
                double[] sampledMutationRates = new double[numberOfFeatures];

                for (int k = 0; k < numberOfFeatures; k++) {
                    sampledMutationRates[k] = solutions[random.nextInt(mu)].mutationRates[k];
                }

                // top individual generates offspring
                newSolutions[i] = solutions[j].generateOffspring(sampledMutationRates);

                if (j++ > mu) {
                    // cycle top mu solutions
                    j = 0;
                }
            }

            return new SolutionSet(newSolutions);
        }

        void evaluate() {
            Arrays.stream(solutions).parallel().forEach(Solution::evaluate);
        }

        private void swap(int i, int j) {
            Solution tmp = solutions[i];
            solutions[i] = solutions[j];
            solutions[j] = tmp;
        }

        public Solution[] getSolutions() {
            return solutions;
        }

        public Solution getBestSolution() {
            return solutions[0];
        }
    }

    public class Solution {
        private final double[] features, mutationRates;
        private double fitness, penalty;

        Solution() {
            this.features = new double[numberOfFeatures];
            this.mutationRates = new double[numberOfFeatures];
        }

        Solution(double[] features, double[] mutationRates) {
            this.features = features;
            this.mutationRates = mutationRates;
        }

        public double[] getFeatures() {
            return features;
        }

        public double getFitness() {
            return fitness;
        }

        public double getPenalty() {
            return penalty;
        }

        Solution generateOffspring(double[] sampledMutationRates) {
            Solution offspring = new Solution();

            double globalLearningRate = tauDash * random.nextGaussian();

            for (int i = 0; i < numberOfFeatures; i++) {
                // global intermediate recombination for mutation rates
                double newMutationRate = (mutationRates[i] + sampledMutationRates[i]) / 2 *
                        // lognormal update for mutation rates
                        Math.exp(globalLearningRate + tau * random.nextGaussian());
                offspring.mutationRates[i] = newMutationRate;
                // generate offspring features
                offspring.features[i] = features[i]; // in case fail to bound

                double newFeatureValue;
                for (int j = 0; j < boundingTries; j++) {
                    // try generating new feature value and check whether it is in bounds
                    newFeatureValue = features[i] + newMutationRate * random.nextGaussian();
                    if (objective.inBounds(newFeatureValue, i)) {
                        offspring.features[i] = newFeatureValue;
                        break;
                    }
                }

            }

            return offspring;
        }

        void evaluate() {
            Objective.Result result = objective.evaluate(features);
            fitness = result.getValue();
            penalty = result.getPenalty();
        }
    }
}
