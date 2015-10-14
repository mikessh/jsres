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

import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

import static java.lang.Math.*;

public class BenchmarkTest {
    @Test
    public void bananaTest() {
        testSres(100, new TestProblem(new double[]{-10, -10}, new double[]{10, 10}, false, 0) {
            @Override
            public Result evaluate(double[] features) {
                double a = (1 - features[0]), b = 10 * (features[1] - features[0] * features[0]);
                return new Result(a * a + b * b);
            }
        });
    }

    @Test
    public void moreBananaTest() {
        int N = 10;
        double[] lowerBounds = new double[N],
                upperBounds = new double[N];

        Arrays.setAll(lowerBounds, value -> -10);
        Arrays.setAll(upperBounds, value -> 10);

        testSres(5000, new TestProblem(lowerBounds, upperBounds, false, 0) {
            @Override
            public Result evaluate(double[] features) {
                double sum = 0;
                for (int i = 0; i < features.length - 1; i++) {
                    double a = (1 - features[i]), b = 10 * (features[i + 1] - features[i] * features[i]);
                    sum += a * a;
                    sum += b * b;
                }
                return new Result(sum / N);
            }
        });
    }

    @Test
    public void g08Test() {
        testSres(200, new TestProblem(new double[]{0, 0}, new double[]{10, 10}, true, 0.095825) {
            @Override
            public Result evaluate(double[] features) {
                double f = pow(sin(2 * PI * features[0]), 3) * sin(2 * PI * features[1]) /
                        (features[0] * features[0] * features[0] * (features[0] + features[1])),
                        g1 = features[0] * features[0] - features[1] + 1,
                        g2 = 1 - features[0] + (features[1] - 4) * (features[1] - 4);
                return new Result(f, new double[]{g1, g2});
            }
        });
    }

    public void testSres(int numberOfGenerations, TestProblem testProblem) {
        SRES sres = new SRES(testProblem);
        SRES.Solution solution = sres.run(numberOfGenerations).getBestSolution();
        System.out.println(solution.toString());
        Assert.assertTrue("Problem is solved up to required precision and constraint violation rules.",
                testProblem.isSolved(solution));
    }

    public static abstract class TestProblem extends Objective {
        private final double valueAtSolution;
        private final double absolutePrecision, relativePrecision;
        private final int allowedViolatedConstraintsCount;

        public TestProblem(double[] featureUpperBounds, double[] featureLowerBounds, boolean maximizationProblem,
                           double valueAtSolution) {
            this(featureUpperBounds, featureLowerBounds, maximizationProblem,
                    valueAtSolution, 0);
        }

        public TestProblem(double[] featureUpperBounds, double[] featureLowerBounds, boolean maximizationProblem,
                           double valueAtSolution, int allowedViolatedConstraintsCount) {
            this(featureUpperBounds, featureLowerBounds, maximizationProblem,
                    valueAtSolution, allowedViolatedConstraintsCount, 1e-6, 1e-3);
        }

        public TestProblem(double[] featureUpperBounds, double[] featureLowerBounds, boolean maximizationProblem,
                           double valueAtSolution, int allowedViolatedConstraintsCount,
                           double absolutePrecision, double relativePrecision) {
            super(featureUpperBounds, featureLowerBounds, maximizationProblem);
            this.valueAtSolution = valueAtSolution;
            this.absolutePrecision = absolutePrecision;
            this.relativePrecision = relativePrecision;
            this.allowedViolatedConstraintsCount = allowedViolatedConstraintsCount;
        }

        public boolean isSolved(SRES.Solution solution) {
            Objective.Result result = solution.getObjectiveResult();

            int violatedConstraints = 0;

            for (int i = 0; i < result.getConstraintValues().length; i++) {
                if (result.getConstraintValues()[i] > 0) {
                    violatedConstraints++;
                }
            }

            double absoluteError = abs(valueAtSolution - result.getValue()),
                    relativeError = absoluteError / abs(valueAtSolution);

            return violatedConstraints <= allowedViolatedConstraintsCount &&
                    absoluteError <= absolutePrecision &&
                    (abs(valueAtSolution) < relativePrecision || relativeError <= relativePrecision);

        }
    }

}
