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

/**
 * A class that implements Stocastic Ranking Evolution Strategy (SRES),
 * an evolutionary algorithm for constrained optimization of real multivariate objective functions.
 * User should provide an objective function instance inherited from abstract {@link com.antigenomics.jsres.Objective} class.
 * The algorithm is executed via {@link #run} method. Objective function evaluation is optimized by implementing
 * it in parallel for all evolved solutions in {@link com.antigenomics.jsres.SRES.SolutionSet}.
 *
 * @see <a href="https://notendur.hi.is/~tpr/software/sres/Tec311r.pdf">Runarsson TP and Yao X. Stochastic Ranking for Constrained Evolutionary Optimization. IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, 2000.</a>
 */
public final class SRES {
    static final Random random = new Random(51102);
    private static final int defaultNumberOfGenerations = 1000, boundingTries = 10;
    private final Objective objective;
    private final double tau, tauDash, rankingPenalizationFactor;
    private final int lambda, mu, numberOfFeatures, numberOfSweeps;

    private final boolean verbose;

    /**
     * Creates an instance of SRES algorithm for a given objective function.
     * Default algorithm parameters from Runarsson TP and Yao X work are used.
     *
     * @param objective objective function instance
     */
    public SRES(Objective objective) {
        this(objective, true);
    }

    /**
     * Creates an instance of SRES algorithm for a given objective function.
     * Default algorithm parameters from Runarsson TP and Yao X work are used.
     *
     * @param objective objective function instance
     * @param verbose   if set to true, will print the number of generations passed and other statistics to {@code stderr}
     */
    public SRES(Objective objective, boolean verbose) {
        this(objective, verbose, 200, 30, 1.0, 200, 0.45);
    }

    /**
     * Creates an instance of SRES algorithm for a given objective function.
     *
     * @param objective                 objective function instance
     * @param verbose                   if set to true, will print the number of generations passed and other statistics to {@code stderr}
     * @param lambda                    solution set size
     * @param mu                        number of top-ranking solutions selected to produce new solution set at each generation
     * @param expectedConvergenceRate   expected convergence rate
     * @param numberOfSweeps            number of times stochastic ranking bubble-sort is applied to solution set
     * @param rankingPenalizationFactor constraint breaking penalization factor, should be in {@code [0, 1]};
     *                                  no penalization is performed if set to 1, all solutions will be penalized if set to 0
     */
    public SRES(Objective objective, boolean verbose,
                int lambda, int mu,
                double expectedConvergenceRate,
                int numberOfSweeps, double rankingPenalizationFactor) {
        this.objective = objective;
        this.verbose = verbose;
        this.numberOfFeatures = objective.getNumberOfFeatures();
        this.lambda = lambda;
        this.mu = mu;
        this.tau = expectedConvergenceRate / Math.sqrt(2 * Math.sqrt(numberOfFeatures));
        this.tauDash = expectedConvergenceRate / Math.sqrt(2 * numberOfFeatures);
        this.numberOfSweeps = numberOfSweeps;
        this.rankingPenalizationFactor = rankingPenalizationFactor;
    }

    /**
     * Runs the SRES algorithm for {@value #defaultNumberOfGenerations} generations.
     *
     * @return a sorted (from best to worst) set of solutions obtained after running the algorithm
     */
    public SolutionSet run() {
        return run(1000);
    }

    /**
     * Runs the SRES algorithm for {@code numberOfGenerations} generations.
     *
     * @param numberOfGenerations number of generations to run
     * @return a sorted (from best to worst) set of solutions obtained after running the algorithm
     */
    public SolutionSet run(int numberOfGenerations) {
        SolutionSet solutionSet = new SolutionSet();

        for (int i = 0; i < numberOfGenerations; i++) {
            solutionSet.evaluate();
            solutionSet.sort();

            if (verbose && i % 100 == 0) {
                Solution bestSolution = solutionSet.getBestSolution();
                System.err.println("SRES ran for " + i + " generations, best solution fitness is " +
                        bestSolution.getFitness() + ", penalty is " + bestSolution.getPenalty());
            }

            solutionSet = solutionSet.evolve();
        }

        solutionSet.evaluate();
        solutionSet.sort();

        if (verbose) {
            Solution bestSolution = solutionSet.getBestSolution();
            System.err.println("SRES finished, best solution fitness is " +
                    bestSolution.getFitness() + ", penalty is " + bestSolution.getPenalty());
        }

        return solutionSet;
    }

    /**
     * A set of solutions for the objective function and constraints.
     */
    public class SolutionSet {
        private final Solution[] solutions;

        /**
         * Generates initial solution set at random.
         */
        SolutionSet() {
            this.solutions = new Solution[lambda];
            for (int i = 0; i < lambda; i++) {
                solutions[i] = new Solution(objective.generateFeatureVector(), objective.getMutationRates());
            }
        }

        /**
         * Constructs a new instance of solution set.
         *
         * @param solutions solutions array
         */
        SolutionSet(Solution[] solutions) {
            this.solutions = solutions;
        }

        /**
         * Sorts solutions using stochastic ranking bubble sort.
         */
        void sort() {
            for (int i = 0; i < numberOfSweeps; i++) {
                boolean swapped = false;
                for (int j = 0; j < lambda - 1; j++) {
                    double p = random.nextDouble(),
                            p1 = solutions[j].getPenalty(),
                            p2 = solutions[j + 1].getPenalty();
                    if (p < rankingPenalizationFactor || (p1 == 0 && p2 == 0)) {
                        // no solution break constraints or penalization is not applied
                        if (solutions[j].getFitness() < solutions[j + 1].getFitness()) {
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

        /**
         * Evolve the population by selecting top {@link #mu} solutions and generating
         * new population of {@link #lambda} solutions using recombination of mutation rates
         * and random weighted mutation of features.
         *
         * @return evolved solution set, no evaluation of objective function is constraints is performed on this step
         */
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

        /**
         * Evaluates objective function and constraints for the solution set in parallel.
         */
        void evaluate() {
            Arrays.stream(solutions).parallel().forEach(Solution::evaluate);
        }

        /**
         * Swaps two solutions.
         *
         * @param i index of fist solution
         * @param j index of second solution
         */
        private void swap(int i, int j) {
            Solution tmp = solutions[i];
            solutions[i] = solutions[j];
            solutions[j] = tmp;
        }

        /**
         * Gets the solutions array.
         *
         * @return array of solutions in this set
         */
        public Solution[] getSolutions() {
            return solutions;
        }

        /**
         * Gets the best solution from solution set.
         *
         * @return best solution
         */
        public Solution getBestSolution() {
            return solutions[0];
        }
    }

    /**
     * A class representing solution to a problem that consists of objective function and constraints.
     */
    public class Solution {
        private final double[] features, mutationRates;
        private Objective.Result objectiveResult = null;

        /**
         * Creates a solution, features and mutation rates are not set.
         */
        Solution() {
            this.features = new double[numberOfFeatures];
            this.mutationRates = new double[numberOfFeatures];
        }

        /**
         * Creates a solution with specified features and mutation rates.
         * Used to instantiate a descendants of solutions selected by stochastic ranking.
         *
         * @param features      feature array
         * @param mutationRates array of mutation rates
         */
        Solution(double[] features, double[] mutationRates) {
            this.features = features;
            this.mutationRates = mutationRates;
        }

        /**
         * Gets the features (objective function parameters) for this solution.
         *
         * @return objective function parameters
         */
        public double[] getFeatures() {
            return features;
        }

        /**
         * Gets the result of evaluation of objective function and constraints for this solution.
         *
         * @return objective function and constraints evaluation results
         */
        public Objective.Result getObjectiveResult() {
            return objectiveResult;
        }

        /**
         * Gets the fitness of this solution.
         * @return {@code value} of objective function for maximization problem, {@code -value} for minimization problem
         * @see com.antigenomics.jsres.Objective.Result#getValue
         */
        double getFitness() {
            return objectiveResult.getObjective().isMaximizationProblem() ?
                    objectiveResult.getValue() : -objectiveResult.getValue();
        }

        /**
         * Gets the penalty value of this solution. 
         * @return penalty value, sum of squares of constraint values for constraints that are violated
         * @see com.antigenomics.jsres.Objective.Result#getPenalty 
         */
        double getPenalty() {
            return objectiveResult.getPenalty();
        }

        /**
         * Generates an offspring for this solution. Mutation rates are recombined and 
         * applied to randomly mutate the feature vector.
         * @param sampledMutationRates random sample of mutation rates from top {@link #mu} solutions for recombination
         * @return a recombined and mutated offspring
         */
        Solution generateOffspring(double[] sampledMutationRates) {
            Solution offspring = new Solution();

            double globalLearningRate = tauDash * random.nextGaussian();

            double[] mutationRateConstraints = objective.getMutationRates();

            for (int i = 0; i < numberOfFeatures; i++) {
                // global intermediate recombination for mutation rates
                double newMutationRate = (mutationRates[i] + sampledMutationRates[i]) / 2 *
                        // lognormal update for mutation rates
                        Math.exp(globalLearningRate + tau * random.nextGaussian());
                // important: prevent mutation rates from growing to infinity
                offspring.mutationRates[i] = Math.min(newMutationRate, mutationRateConstraints[i]);
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

        /**
         * Evaluates given solution using objective function and constraints specified in parent {@code SRES} algorithm. 
         */
        void evaluate() {
            objectiveResult = objective.evaluate(features);
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public String toString() {
            return "SOLUTION " + Arrays.toString(features) + ",fitness is " + getFitness() + ", penalty is " +
                    getPenalty() + "\n" + (objectiveResult == null ? "OBJECTIVE NOT EVALUATED" : objectiveResult.toString());
        }
    }
}
