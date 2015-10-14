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

public class BenchmarkTest {
    @Test
    public void bananaTest() {
        testSres(new TestProblem(new double[]{-10, -10},
                new double[]{10, 10}, false, 0) {
            @Override
            public Result evaluate(double[] features) {
                double a = (1 - features[0]), b = 10 * (features[1] - features[0] * features[0]);
                return new Result(a * a + b * b);
            }
        });
    }

    public void testSres(TestProblem testProblem) {
        SRES sres = new SRES(testProblem);
        SRES.Solution solution = sres.run().getBestSolution();
        System.out.println(solution.toString());
        Assert.assertTrue("Problem is solved up to required precision and constraint violation rules.",
                testProblem.isSolved(solution));
    }

    public static abstract class TestProblem extends Objective {
        private static final double EPS = 1e-16;
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

            double absoluteError = Math.abs(valueAtSolution - result.getValue()),
                    relativeError = absoluteError / (Math.abs(valueAtSolution) + EPS);

            return violatedConstraints <= allowedViolatedConstraintsCount &&
                    absoluteError <= absolutePrecision && relativeError <= relativePrecision;

        }
    }

}
