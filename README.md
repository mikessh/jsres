### Description

A java-based implementation of Stochastic Ranking Evolutionary Strategy (SRES). 
In a nutshell, this is a heuristic algorithm that finds a set of parameters ``Xi`` that minimize/maximize a given function's value (objective) ``F(Xi)`` and does not violate the constraints specified in 
for of ``Gk(Xi)<=0`` inequalities.
For more details, please see the [Runarsson TP and Yao X paper](https://notendur.hi.is/~tpr/software/sres/Tec311r.pdf) which has served as a basis for implementation of this algorithm. In my experience, this is one of the best algorithms 
for solving real multivariate objective function constrained optimization problems. Note that this implementation is optimized by calculating the objective and constraints in parallel over the population of evolving solutions.

### Installation

[Java 1.8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) is required for this library. In order to use, compile using 
[Maven](https://maven.apache.org/) and add the following dependency to your ``pom.xml`` file:

```xml
<dependency>
	<groupId>com.antigenomics</groupId>
    <artifactId>jsres</artifactId>
    <version>1.0.0</version>            
</dependency>
```

### Usage example

The following code demonstrates how to specify the objective function and constraints:

```java
import static java.lang.Math.*;

...

Objective myObjective = new Objective(
        new double[]{0, 0},   // lower bounds of the search space
        new double[]{10, 10}, // upper bounds of the search space
        true // this is a maximization problem
    ) {
        @Override
        public Result evaluate(double[] features) {
        // objective, first and second constraint values
        double f = pow(sin(2 * PI * features[0]), 3) * sin(2 * PI * features[1]) /
                   (features[0] * features[0] * features[0] * (features[0] + features[1])),
               g1 = features[0] * features[0] - features[1] + 1,
               g2 = 1 - features[0] + (features[1] - 4) * (features[1] - 4);
        return new Result(f, new double[]{g1, g2});
    }
}
```

The following code can be used to execute the SRES algorithm:

```java
int numberOfGenerations = 1000; // run for 1000 generations
SRES sres = new SRES(myObjective); // set objective
SRES.SolutionSet solutionSet = sres.run(numberOfGenerations); // holds solutions after last SRES generation
```

Accessing the results is performed as follows:

```java
SRES.Solution bestSolution = solutionSet.getBestSolution();
System.out.println(solution.toString()); // print the best solution
solution.getFeatures(); // optimal objective parameters
solution.getObjectiveResult().getValue(); // objective function value at optimum
solution.getObjectiveResult().getConstraintValues(); // values of constraints at optimum
```