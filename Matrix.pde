// So many matrix solvers out there
// Not many of them support PCG, and we don't know who performs better solving those systems
// so I'm setting up a wrapper class for matrix
// and hopefully that will make it easy for us to switch between libraries in the future.
import org.apache.commons.math3.linear.*;

/** 
 * Sparse Matrix
 */
 
class SparseMatrix {
  protected final AbstractRealMatrix matrix;
  
  /**
   * Create a new sparse matrix
   *
   * @param number of columns
   * @param number of rows
   */
  public SparseMatrix(int columns, int rows) {
    matrix = new Array2DRowRealMatrix(rows, columns);
  }
  
  public void setEntry(int column, int row, double value) {
    matrix.setEntry(row, column, value); // yes, their rows come before columns
  }
  
  public double getEntry(int column, int row) {
    return matrix.getEntry(row, column); 
  }
}

/**
 * Vector
 */
class Vector {
  protected final RealVector vector;
  
  public Vector(int rows) {
    vector = new ArrayRealVector(rows);
  }
  
  protected Vector(RealVector vector) {
    this.vector = vector; 
  }
  
  public void setEntry(int row, double value) {
    vector.setEntry(row, value);
  }
  
  public double getEntry(int row) {
    return vector.getEntry(row); 
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
//  private final int MAX_PCG_ITERATIONS = 100;
//  private final double PCG_DELTA = 2;
//  private final boolean PCG_CHECK_POSITIVE_DEFINITENESS = true;
//  private final ConjugateGradient conjugateGradientSolver = 
//                      new ConjugateGradient(MAX_PCG_ITERATIONS, // maximum iterations
//                                            PCG_DELTA, 
//                                            PCG_CHECK_POSITIVE_DEFINITENESS);
  
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
  public Vector pcgSolve(SparseMatrix A, Vector b) {
//    Vector initialGuess = new Vector(b.size());
//    for (int i = 0; i < initialGuess.size(); i++)
//      initialGuess.setEntry(i, 0); // not sure what to put here, honestly 
//    ConjugateGradient conjugateGradientSolver = new ConjugateGradient(1000, 100, false);
//    
//    RealVector x = conjugateGradientSolver.solve(A.matrix, b.vector);
//    println(x.getEntry(0));
//    return new Vector(x);

    DecompositionSolver solver = new LUDecomposition(A.matrix).getSolver();
    RealVector x = solver.solve(b.vector);
    println(x.getEntry(0));
    return new Vector(x);
  }
}
