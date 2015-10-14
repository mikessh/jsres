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

public abstract class Objective {
    private final int numberOfFeatures;
    private final double[] featureLowerBounds, featureUpperBounds;
    private final double[] mutationRates;
    private final boolean maximizationProblem;

    protected Objective(double[] featureLowerBounds, double[] featureUpperBounds,
                        boolean maximizationProblem) {
        this.numberOfFeatures = featureLowerBounds.length;
        if (featureUpperBounds.length != numberOfFeatures) {
            throw new IllegalArgumentException("Lengths of upper and lower bounds should match.");
        }
        this.featureLowerBounds = featureLowerBounds;
        this.featureUpperBounds = featureUpperBounds;
        for (int i = 0; i < numberOfFeatures; i++) {
            if (featureUpperBounds[i] < featureLowerBounds[i]) {
                throw new IllegalArgumentException("Feature upper bound is smaller than lower bound for index " + i + ".");
            }
        }
        this.maximizationProblem = maximizationProblem;

        this.mutationRates = new double[numberOfFeatures];
        double factor = Math.sqrt(numberOfFeatures);

        for (int i = 0; i < numberOfFeatures; i++) {
            mutationRates[i] = (featureUpperBounds[i] - featureLowerBounds[i]) / factor;
        }
    }

    public double[] generateFeatureVector() {
        double[] features = new double[numberOfFeatures];

        for (int i = 0; i < numberOfFeatures; i++) {
            double x = SRES.random.nextDouble();
            features[i] = featureLowerBounds[i] + x * (featureUpperBounds[i] - featureLowerBounds[i]);
        }

        return features;
    }

    public final double[] getMutationRates() {
        return mutationRates;
    }

    public final boolean inBounds(double feature, int index) {
        return featureLowerBounds[index] <= feature && feature <= featureUpperBounds[index];
    }

    public abstract Result evaluate(double[] features);

    public final int getNumberOfFeatures() {
        return numberOfFeatures;
    }

    public boolean isMaximizationProblem() {
        return maximizationProblem;
    }

    public class Result {
        private final double value, penalty;
        private final double[] constraintValues;

        public Result(double value) {
            this.value = value;
            this.constraintValues = new double[0];
            this.penalty = 0;
        }

        public Result(double value, double[] constraintValues) {
            this.value = value;
            this.constraintValues = constraintValues;
            double penalty = 0;
            for (double p : constraintValues) {
                double pp = Math.max(0, p);
                penalty += pp * pp;
            }
            this.penalty = penalty;
        }

        public double[] getConstraintValues() {
            return constraintValues;
        }

        public double getValue() {
            return value;
        }

        double getPenalty() {
            return penalty;
        }

        public Objective getObjective() {
            return Objective.this;
        }

        @Override
        public String toString() {
            return "RESULT objective function value is " + value + ", constraint values are " +
                    Arrays.toString(constraintValues);
        }
    }
}
