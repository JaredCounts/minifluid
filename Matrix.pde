// So many matrix solvers out there
// Not many of them support PCG, and we don't know who performs better solving those systems
// so I'm setting up a wrapper class for matrix
// and hopefully that will make it easy for us to switch between libraries in the future.
import org.apache.commons.math3.linear.*;

/** 
 * Matrix
 *
 * Immutable object representing a matrix.
 */
 
class Matrix {
  protected final RealLinearOperator matrix;
  
  /**
   * Create a new matrix from a 2D double array
   *
   * @param doubleMatrix
   *     Matrix of doubles where doubleMatri[j][i] gives us the cell at row j, column i.
   *     Must have at least one element.
   */
  public Matrix(double[][] doubleMatrix) {
    matrix = new Array2DRowRealMatrix(doubleMatrix);
  }
  
  protected Matrix(RealLinearOperator matrix) {
    this.matrix = matrix; 
  }
  
}

/**
 * Vector
 *
 * Immutable object representing a column vector.
 */
class Vector {
  protected final RealVector vector;
  
  public Vector(double[] doubleVector) {
    vector = new ArrayRealVector(doubleVector);
  }
  
  protected Vector(RealVector vector) {
    this.vector = vector; 
  }
  
  public int size() {
    return vector.getDimension();
  }
}

/**
 * Linear Solver
 *
 * Provides helper functions for solving linear systems of the form Ax = b
 * Where A is a 2D matrix, x is the vector of unknowns, and b is a vector.
 * 
 */
class LinearSolver {
  private final int MAX_PCG_ITERATIONS = 100;
  private final double PCG_DELTA = 0.01;
  private final boolean PCG_CHECK_POSITIVE_DEFINITENESS = true;
  private final ConjugateGradient conjugateGradientSolver = 
                      new ConjugateGradient(MAX_PCG_ITERATIONS, // maximum iterations
                                            PCG_DELTA, 
                                            PCG_CHECK_POSITIVE_DEFINITENESS);
  
  /**
   * PCG Solver
   *
   * Solves Ax = b using the preconditioned conjugate gradient method.
   *
   * @param Matrix A from Ax = b. Must be symmetric and positive definite.
   * 
   * @param Vector b of known values from Ax = b.
   *
   */
  public Vector pcgSolve(Matrix A, Vector b) {
    double[] initialGuess = new double[b.size()];
    for (int i = 0; i < initialGuess.length; i++)
      initialGuess[i] = Math.random(); // not sure what to put here, honestly 
    
    Vector x0 = new Vector(initialGuess);
    
    RealVector x = conjugateGradientSolver.solveInPlace(A.matrix,
                                                        null,
                                                        b.vector,
                                                        x0.vector);
    return new Vector(x);
  }
}
