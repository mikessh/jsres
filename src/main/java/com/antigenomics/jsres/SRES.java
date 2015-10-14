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

    public Population run(int numberOfGenerations) {
        Population population = new Population();

        for (int i = 0; i < numberOfGenerations; i++) {
            population.evaluate();
            population.sort();
            population = population.evolve();
        }

        population.evaluate();

        return population;
    }

    public class Population {
        private final Individual[] individuals;

        public Population() {
            this.individuals = new Individual[lambda];
            for (int i = 0; i < lambda; i++) {
                individuals[i] = new Individual(objective.generateFeatureVector(), objective.getMutationRates());
            }
        }

        public Population(Individual[] individuals) {
            this.individuals = individuals;
        }

        public void sort() {
            for (int i = 0; i < numberOfSweeps; i++) {
                boolean swapped = false;
                for (int j = 0; j < lambda - 1; j++) {
                    double p = random.nextDouble(),
                            p1 = individuals[j].getPenalty(),
                            p2 = individuals[j + 1].getPenalty();
                    if (p < rankingPenalizationFactor || (p1 == 0 && p2 == 0)) {
                        if (individuals[j].getFitness() > individuals[j + 1].getFitness()) {
                            swap(j, j + 1);
                            swapped = true;
                        }
                    } else if (p1 > p2) {
                        swap(j, j + 1);
                        swapped = true;
                    }
                }
                if (!swapped) {
                    break;
                }
            }
        }

        public Population evolve() {
            Individual[] newIndividuals = new Individual[lambda];

            int j = 0;
            for (int i = 0; i < lambda; i++) {
                // mutation rates are sampled from top individuals
                double[] sampledMutationRates = new double[numberOfFeatures];

                for (int k = 0; k < numberOfFeatures; k++) {
                    sampledMutationRates[k] = individuals[random.nextInt(mu)].getMutationRate(k);
                }

                // top individual generates offspring
                newIndividuals[i] = individuals[j].generateOffspring(sampledMutationRates);

                if (j++ > mu) {
                    // cycle top mu individuals
                    j = 0;
                }
            }

            return new Population(newIndividuals);
        }

        public void evaluate() {
            Arrays.stream(individuals).parallel().forEach(Individual::evaluate);
        }

        private void swap(int i, int j) {
            Individual tmp = individuals[i];
            individuals[i] = individuals[j];
            individuals[j] = tmp;
        }
    }

    public class Individual {
        private final double[] features, mutationRates;
        private double fitness, penalty;

        public Individual() {
            this.features = new double[numberOfFeatures];
            this.mutationRates = new double[numberOfFeatures];
        }

        public Individual(double[] features, double[] mutationRates) {
            this.features = features;
            this.mutationRates = mutationRates;
        }

        public double getFeature(int index) {
            return features[index];
        }

        public double getMutationRate(int index) {
            return mutationRates[index];
        }

        public double getFitness() {
            return fitness;
        }

        public double getPenalty() {
            return penalty;
        }

        public Individual generateOffspring(double[] sampledMutationRates) {
            Individual offspring = new Individual();

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

        public void evaluate() {
            Objective.Result result = objective.compute(features);
            fitness = result.getValue();
            penalty = result.getPenalty();
        }
    }
}
