/** 
 * Fluid Grid
 * Mutable object.
 * MAC staggered grid where velocities are stored at cell edges and pressure in cell centers.
 * Also contains helper methods and a solve method for moving the grid values forward in time.
 */
class FluidGrid {
  // cells will be square
  private final float cellWidth;

  // grid of cells
  // to access, use cell cells[i][j], where cell is at the i'th column and j'th row.
  // cells[0][0]'s top left corner is at (0,0).
  private final FluidGridCell[][] cells;
  
  // when updating velocities, we will use these buffers to store computed values
  // this is so updated cells don't corrupt computation for cells not yet updated
  // keeping one buffer to be used everywhere helps with caching and memory management
  // that is, our code can just reuse this memory instead of allocating new chunks for each step.
  float[][] velocityXBuffer;
  float[][] velocityYBuffer;

  // the solver we'll be using to solve for incompressability
  private LinearSolver solver;

  /** 
   * Constructor for fluid grid
   * Initializes all of the cells in the grid.
   * For best results, choose a cellWidth that's a common factor of regionWidth and regionHeight.
   * @param cellWidth height and width of each square cell
   *                  must be > 0
   * @param regionWidth width of the region to solve fluid in
   * @param regionHeight height of region
   */
  public FluidGrid(float cellWidth, float regionWidth, float regionHeight) {
    solver = new LinearSolver();

    jAssert("FluidGrid: cellWidth <= 0", cellWidth > 0);
    jAssert("FluidGrid: regionWidth < 0", regionWidth > 0);
    jAssert("FluidGrid: regionHeight < 0", regionHeight > 0);

    this.cellWidth = cellWidth;

    int columnCount = ceil(regionWidth / cellWidth);
    int rowCount = ceil(regionHeight / cellWidth);

    cells = new FluidGridCell[columnCount][rowCount];
    
    velocityXBuffer = new float[cells.length+1][cells[0].length+1];
    velocityYBuffer = new float[cells.length+1][cells[0].length+1];
    
    for (int i = 0; i < columnCount; i++) {
      for (int j = 0; j < rowCount; j++) {
        FluidGridCell cell = new FluidGridCell();
        cells[i][j] = cell;
        
        cell.hasLiquid = true; //random(0,1) > 0.5;

        cell.position = new PVector(i * cellWidth + cellWidth/2.f, j * cellWidth + cellWidth/2.f);
        cell.topEdgePosition = PVector.add(cell.position, new PVector(0, -cellWidth/2.f));
        cell.rightEdgePosition = PVector.add(cell.position, new PVector(cellWidth/2.f,0));
        cell.bottomEdgePosition = PVector.add(cell.position, new PVector(0,cellWidth/2.f));
        cell.leftEdgePosition = PVector.add(cell.position, new PVector(-cellWidth/2.f,0));

        // initialize left/top edge velocities
        Float velocityX = new Float(random(-10, 10));
        cell.velocityXLeft = velocityX;
        if (i != 0) // adjacent cells will point to the same velocity objects
          cells[i-1][j].velocityXRight = velocityX;

        Float velocityY = new Float(random(-10, 10));
        cell.velocityYTop = velocityY;
        if (j != 0)
          cells[i][j-1].velocityYBottom = velocityY;

        // also take care of right/bottom edge velocities for far right/bottom cells
        if (i == columnCount-1)
          cell.velocityXRight = new Float(0);

        if (j == rowCount-1)
          cell.velocityYBottom = new Float(0);
      }
    }
  }

  /**
   * Get the fluid grid's region width
   */
  public float getRegionWidth() {
    return cellWidth * cells.length;
  }
  /**
   * Get the fluid grid's region height
   */
  public float getRegionHeight() {
    return cellWidth * cells[0].length;
  }


  public void draw() {
    for (FluidGridCell[] cellColumn : cells) {
      for (FluidGridCell cell : cellColumn)
        cell.draw();
    }
  }
  
  /**
   * This will set the edge velocity values of cells to 
   * their updated values found in their corresponding positions velocityXBuffer and velocityYBuffer
   * this will also reset the buffers for later usage.
   */
  void loadVelocityBuffersIntoCells() {
    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[i].length; j++) {
        FluidGridCell cell = cells[i][j];

        cell.velocityXLeft = velocityXBuffer[i][j];
        cell.velocityYTop  = velocityYBuffer[i][j];
        
        velocityXBuffer[i][j] = 0;
        velocityYBuffer[i][j] = 0;

        if (i == cells.length-1) {
          cell.velocityXRight  = velocityXBuffer[i+1][j];
          velocityXBuffer[i+1][j] = 0;  
        }
        if (j == cells[i].length-1) {
          cell.velocityYBottom = velocityYBuffer[i][j+1];
          velocityYBuffer[i][j+1] = 0;
        }
      }
    }
  }
  
  void applyExternalForces(float timestep) {
    // just integrate gravity into cell edge velocities
    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[i].length; j++) {
        FluidGridCell cell = cells[i][j];

        // since edges are shared between cells, we only solve for left and top
        cell.velocityXLeft += GRAVITY.x * timestep; // a bit redundant, but lets not make assumptions
        cell.velocityYTop  += GRAVITY.y * timestep;

        // except for the bottom-most and right-most cells
        if (i == cells.length-1)
          cell.velocityXRight  += GRAVITY.x * timestep;
        if (j == cells[i].length-1)
          cell.velocityYBottom += GRAVITY.y * timestep;
      }
    }
  }
  
  void applyConvection(float timestep) {
    // using semi-lagrangian method
    //   Stam, "Stable Fluids"
    //   http://www.autodeskresearch.com/pdf/ns.pdf
    // for any given cell, make an imaginary particle in the middle
    // then we trace this particle back in time by timestep
    //   do this using Runge-Kutta or some other higher-order method
    // then take the "Averaged" velocity at that point and move it to the current cell

    // however, cells don't *have* a velocity vector at the center
    // so do we do component-only solves for the centers of edges for now?
    // Sure.
    
    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[0].length; j++) {
        FluidGridCell cell = cells[i][j];

        // solve for left edge
        PVector leftEdgePosition = cell.leftEdgePosition;
        PVector leftEdgeTracedPosition = traceFrom(leftEdgePosition, -timestep);
        PVector leftEdgeTracedVelocity = getVelocityAt(leftEdgeTracedPosition);
        velocityXBuffer[i][j] = leftEdgeTracedVelocity.x;

        // solve for upper edge
        PVector topEdgePosition = cell.topEdgePosition;
        PVector topEdgeTracedPosition = traceFrom(topEdgePosition, -timestep);
        PVector topEdgeTracedVelocity = getVelocityAt(topEdgeTracedPosition);
        velocityYBuffer[i][j] = topEdgeTracedVelocity.y;

        // lets not forget the last cells' bottom and right velocities
        if (i == cells.length-1 && j < cells[0].length) {
          PVector rightEdgePosition = cell.rightEdgePosition;
          PVector rightEdgeTracedPosition = traceFrom(rightEdgePosition, -timestep);
          PVector rightEdgeTracedVelocity = getVelocityAt(rightEdgeTracedPosition);
          velocityXBuffer[i+1][j] = rightEdgeTracedVelocity.x;
        }
        if (j == cells[0].length-1 && i < cells.length) {
          PVector bottomEdgePosition = cell.bottomEdgePosition;
          PVector bottomEdgeTracedPosition = traceFrom(bottomEdgePosition, -timestep);
          PVector bottomEdgeTracedVelocity = getVelocityAt(bottomEdgeTracedPosition);
          velocityYBuffer[i][j+1] = bottomEdgeTracedVelocity.y;
        }
      }
    }

    loadVelocityBuffersIntoCells();
  }
  
  void applyViscosity(float timestep) {
    // using standard central differencing
    //   Foster and Metaxas, "Realistic Animation of Liquids"
    //   http://graphics.stanford.edu/courses/cs468-05-fall/Papers/foster-metaxas-gmip96.pdf

    // for each axis, we add this to velocity:
    // timestep * viscosity * (velocityToLeft - 2*velocity + velocityToRight) / (cellWidth^2)
    // where velocityToLeft and velocityToRight are the velocities one cell over left and right respectively

    // it seems to me that by "averaging" adjacent velocities into the middle by a multiple of viscosity,
    // we're simulating viscous drag. An averaged value will be more resistant to change than any single value.

    // for the same reason as in convection, we'll be saving velocities into a temporary matrix before applying them
    // so that there isn't any new vs old value contamination as we move across each cell
    // we'll be reusing the float 2D arrays newVelocitiesX and newVelocitiesY

    // only need to compute this once
    float viscousFactor = VISCOSITY * timestep / sq(cellWidth);

    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[i].length; j++) {
        FluidGridCell cell = cells[i][j];

        // solve for left edge
        if (i == 0) {
          // there is no left cell
          // so just take the average between the one nearest velocity
          velocityXBuffer[i][j] = cell.velocityXLeft + viscousFactor * (-cell.velocityXLeft + cell.velocityXRight);
        }
        else {
          FluidGridCell leftCell = cells[i-1][j];
          velocityXBuffer[i][j] = cell.velocityXLeft + viscousFactor * (leftCell.velocityXLeft - 2f * cell.velocityXLeft + cell.velocityXRight);
        }

        // solve for upper edge
        if (j == 0) {
          velocityYBuffer[i][j] = cell.velocityYTop + viscousFactor * (-cell.velocityYTop + cell.velocityYBottom);
        }
        else {
          FluidGridCell topCell = cells[i][j-1];
          velocityYBuffer[i][j] = cell.velocityYTop + viscousFactor * (topCell.velocityYTop - 2f * cell.velocityYTop + cell.velocityYBottom);
        }

        // and lets not forget far right and far bottom edge velocities!
        if (i == cells.length-1) {
          velocityXBuffer[i+1][j] = cell.velocityXRight + viscousFactor * (cell.velocityXLeft - cell.velocityXRight);
        }
        if (i == cells.length-1) {
          velocityYBuffer[i][j+1] = cell.velocityYBottom + viscousFactor * (cell.velocityYTop - cell.velocityYBottom);
        }
      }
    }

    loadVelocityBuffersIntoCells();
  }

  void applyIncompressibility(float timestep) {
    // using laplacian + linear system solver
    //   Foster and Fedkiw, "Practical Animation of Liquids"
    //   http://physbam.stanford.edu/~fedkiw/papers/stanford2001-02.pdf
    // Foster and Fedkiw's paper describes how we couple pressure with incompressability via spatial variation
    // 
    // the crucial equation we're solving is the following (per cell):
    // (sum of adjacent pressures) - 4*pressure = (density * cellWidth / timestep) * totalInwardVelocity
    //
    // what we want to do is set up a matrix expression A p = b
    //
    // where p is the vector of unknown pressures to make the velocity field divergence free
    //
    // b is the right hand side of the equation
    //
    // and A is a matrix s.t. a(i,i) = -number of adjacent liquid cells to cell i
    //                        a(i,j) = a(j,i) = 1 for all liquid cells j adjacent to cell i
    // A gives us the coefficients of the pressures on the left hand side
    // and p are the pressure values on the left hand side 
    //
    // cell index in the matrix = (x + width * y + width * height * z) for 3D, or in our case: (x + width * y) for 2D

    // set up liquid matrix and divergence vector
    // which respectively make up the A and b of our A x = b equation.
    // where column or row i corresponds to cellColumn + cellRow * cellColumnCount
    // a(i,i) = negative number of adjacent liquid cells
    // a(i,j) = a(j,i) = 1 if cell j adjacent to cell i is liquid, 0 otherwise
    int matrixWidth = (int)(cells.length * cells[0].length);

    // matrixWidth can be on the order of millions or billions, even
    // yet we only need to put values along the diagonal and off diagonals

    // so it's crucial that we use a sparse matrix,
    // unless we want to use up your computer's entire memory, and take minutes per timestep

    // the liquid matrix gives us coefficients for the left hand side of the equation above
    SparseMatrix liquidMatrix = new SparseMatrix(matrixWidth, matrixWidth);

    Vector divergenceVector = new Vector(matrixWidth);

    // only need to compute this once
    float divergenceFactor = DENSITY * cellWidth / timestep;

    // generating divergenceVector and liquidMatrix
    println("generating divergenceVector and liquidMatrix: " + matrixWidth);
    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[i].length; j++) {
        FluidGridCell cell = cells[i][j];

        int liquidMatrixIndex = i + cells.length * j;
                              
        int adjacentLiquidCellCount = 0;
        if (i != 0) {
          if (cells[i-1][j].hasLiquid)
            adjacentLiquidCellCount++;
        }
        if (i != cells.length-1) {
          if (cells[i+1][j].hasLiquid)
            adjacentLiquidCellCount++;
        }
        if (j != 0) {
          if (cells[i][j-1].hasLiquid)
          adjacentLiquidCellCount++;
        }
        if (j != cells[0].length-1) {
          if (cells[i][j+1].hasLiquid)
            adjacentLiquidCellCount++;
        }
        
        liquidMatrix.setEntry(liquidMatrixIndex, // column
                              liquidMatrixIndex, // row
                              -adjacentLiquidCellCount); // value
                                
        // set off-diagonal values of the pressure matrix
        if (cell.hasLiquid) {
          if (i != 0) {
            if (cells[i-1][j].hasLiquid) {
              int liquidMatrixAdjacentIndex = (i-1) + cells.length * j;
              liquidMatrix.setEntry(liquidMatrixIndex, // column
                                    liquidMatrixAdjacentIndex, // row
                                    1.0); // value
              liquidMatrix.setEntry(liquidMatrixAdjacentIndex, // column
                                    liquidMatrixIndex, // row
                                    1.0); // value
            }
          }
          if (i != cells.length-1) {
            if (cells[i+1][j].hasLiquid) {
              int liquidMatrixAdjacentIndex = (i+1) + cells.length * j;
              liquidMatrix.setEntry(liquidMatrixIndex, // column
                                    liquidMatrixAdjacentIndex, // row
                                    1.0); // value
              liquidMatrix.setEntry(liquidMatrixAdjacentIndex, // column
                                    liquidMatrixIndex, // row
                                    1.0); // value
            }
          }
          if (j != 0) {
            if (cells[i][j-1].hasLiquid) {
              int liquidMatrixAdjacentIndex = i + cells.length * (j-1);
              liquidMatrix.setEntry(liquidMatrixIndex, // column
                                    liquidMatrixAdjacentIndex, // row
                                    1.0); // value
              liquidMatrix.setEntry(liquidMatrixAdjacentIndex, // column
                                    liquidMatrixIndex, // row
                                    1.0); // value
            }
          }
          if (j != cells[0].length-1) {
            if (cells[i][j+1].hasLiquid) {
              int liquidMatrixAdjacentIndex = i + cells.length * (j+1);
              liquidMatrix.setEntry(liquidMatrixIndex, // column
                                    liquidMatrixAdjacentIndex, // row
                                    1.0); // value
              liquidMatrix.setEntry(liquidMatrixAdjacentIndex, // column
                                    liquidMatrixIndex, // row
                                    1.0); // value
            }
          }
        }
        
        // and finally the divergence vector entry
        float totalInwardVelocity = cell.velocityYBottom - cell.velocityYTop + cell.velocityXRight - cell.velocityXLeft;
        float divergence = totalInwardVelocity * divergenceFactor;
        divergenceVector.setEntry(liquidMatrixIndex, (double)divergence);
      }
    }

    // solve for pressures
    // where liquidMatrix * pressures = divergence
    // we use PCG, as suggested by Foster and Fedkiw in their paper
    // since the liquidMatrix is symmetric, sparse, and positive definite
    println("solving for pressures");
    Vector pressures = solver.pcgSolve(liquidMatrix, divergenceVector);

    // update cell pressures 
    println("update cell pressures");
    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[i].length; j++) {
        FluidGridCell cell = cells[i][j];

        int pressuresIndex = i + j * cells.length;

        float pressure = (float)pressures.getEntry(pressuresIndex);

        cell.pressure = pressure;
      }
    }

    // update cell velocities
    println("update cell velocities");
    float pressureToVelocityFactor = timestep / (DENSITY * cellWidth);
    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[i].length; j++) {
        FluidGridCell cell = cells[i][j];

        // update left and top velocities
        if (i != 0)
          cell.velocityXLeft -= pressureToVelocityFactor * (cell.pressure - cells[i-1][j].pressure);
        if (j != 0)
          cell.velocityYTop -= pressureToVelocityFactor * (cell.pressure - cells[i][j-1].pressure);

        // not sure how to account for last bottom and right velocities
      }
    }
  }

  /**
   * Solving
   * Using Foster and Fedkiw's "Practical Animation of Liquids" as reference
   * http://physbam.stanford.edu/~fedkiw/papers/stanford2001-02.pdf
   */
  public void solve() {
    // determine maximum timestep using Courant-Friedrich's-Lewy (CFL) condition
    // as mentioned in Foster and Metaxa's paper,
    // there's a condition that gives us the maximum stable timestep for numerical integration using a course grid
    // 1 > maxAxisVelocity * timestep / cellSize
    // where maxAxisVelocity is the maximum velocity along any edge of a grid cell
    //   thus, timestep < cellSize/maxVelocity
    // for a dynamically sized grid (like a quadtree), it would be
    //   timestep < max(cellSize/cellVelocity)
    // where the right hand side is the maximum of the ratio of cellsize to velocity per each cell.
    float timestep;
    float maxVelocityMagitudeSquared = 0;
    for (FluidGridCell[] cellColumn : cells) {
      for (FluidGridCell cell : cellColumn) {
        // at this point, incompressability should be solved for from the previous solve step
        // therefore the inward velocity of each cell should match the outward velocity
        // so all we need is one x value and one y velocity value from the cell
        // another note, we don't need to do any expensive square-root operations until we find our max value
        float velocityMagnitudeSquared = sq(cell.velocityXLeft) + sq(cell.velocityYTop);
        maxVelocityMagitudeSquared = max(maxVelocityMagitudeSquared, velocityMagnitudeSquared);
      }
    }
    float maxVelocityMagnitude = sqrt(maxVelocityMagitudeSquared);

    // and finally, CFL condition
    timestep = cellWidth / maxVelocityMagnitude;
    jAssert("solve: timestep <= 0", timestep > 0);

    int solveCount;
    if (REAL_TIME) {
      float timeToSolveFor = timeKeeper.getTimeToSolveFor();
      // separate time to solve for into discrete steps
      solveCount = floor(timeToSolveFor / timestep);

      // compute left over time we couldn't account for
      float leftOverTime = timeToSolveFor - solveCount * timestep;

      // store it over for the next solve
      timeKeeper.setLeftOverSolveTime(leftOverTime);
    }
    else {
      solveCount = 1;
      timestep = 0.1;
    }
    
    println("timestep " + timestep + ", solveCount " + solveCount);
    for (int solve = 0; solve < solveCount; solve++) {
      /* -------------- External Forces -------------- */
      println("external forces");
      applyExternalForces(timestep);

      /* -------------- Convection -------------- */
      println("convection");
      applyConvection(timestep);

      /* -------------- Viscosity -------------- */
      println("viscosity");
      applyViscosity(timestep);

      /* -------------- Incompressability and Pressure -------------- */
      println("incompressability");
      applyIncompressibility(timestep);
    }
  }

  /**
   * Returns the resulting position after amountOfTime starting from the given position
   */
  private PVector traceFrom(PVector position, float amountOfTime) {
    // we use 2nd order Runge-Kutta
    // which  means we trace only half way before resampling velocity

    // compute velocity at current location
    PVector velocity = getVelocityAt(position);

    // move along that velocity by half a timestep
    PVector secondPosition = PVector.add(position, PVector.mult(velocity, amountOfTime/2f));
    
    // get velocity to where we've moved back
    PVector secondVelocity = getVelocityAt(secondPosition);

    // move the rest of the way using our updated velocity
    PVector finalPosition = PVector.add(secondPosition, PVector.mult(secondVelocity, amountOfTime/2f));

    // if our boundary conditions are right, finalPosition should still be within the fluid grid
//    jAssert("traceFrom: finalPosition.x < 0", finalPosition.x >= 0);
//    jAssert("traceFrom: finalPosition.x > regionWidth", finalPosition.x < cellWidth * cells.length);
//    jAssert("traceFrom: finalPosition.y < 0", finalPosition.y >= 0);
//    jAssert("traceFrom: finalPosition.y > regionHeight", finalPosition.y < cellWidth * cells[0].length);

    return finalPosition;
  }

  /**
   * Returns a velocity at the given position
   * @param position
   *         must be within grid region
   */
  private PVector getVelocityAt(PVector position) {
    // Since velocity is stored in discrete positions (the cell edges)
    // we must bilinearly interpolate between these positions to get a velocity at other arbitrary positions

    // looking at https://cs.uwaterloo.ca/~c2batty/teaching/general/Fluids2.pdf
    // page 21, we see that the x and y velocities each form their own staggered grids
    // (Hence the name, "MAC staggered grid")
    int column = columnFromPosition(position);
    int row = rowFromPosition(position);
    FluidGridCell containingCell = getCellContaining(position);

    boolean leftHalf = position.x < containingCell.position.x;
    boolean topHalf = position.y < containingCell.position.y;

    FluidGridCell topLeftCell, topRightcell, bottomLeftCell, bottomRightCell;
    
    int leftColumn = max(column - 1, 0);
    int rightColumn = min(column + 1, cells.length - 1);
    int topRow = max(row - 1, 0);
    int bottomRow = min(row + 1, cells[0].length - 1);
    
    if (leftHalf && topHalf) { // top left
      topLeftCell = cells[leftColumn][topRow];
      topRightcell = cells[column][topRow];
      bottomRightCell = cells[column][row];
      bottomLeftCell = cells[leftColumn][row];
    }
    else if (!leftHalf && topHalf) { // top right
      topLeftCell = cells[column][topRow];
      topRightcell = cells[rightColumn][topRow];
      bottomRightCell = cells[rightColumn][row];
      bottomLeftCell = cells[column][row];
    }
    else if (!leftHalf && !topHalf) { // bottom right
      topLeftCell = cells[column][row];
      topRightcell = cells[rightColumn][row];
      bottomRightCell = cells[rightColumn][bottomRow];
      bottomLeftCell = cells[column][bottomRow];
    }
    else if (leftHalf && !topHalf) { // bottom left
      topLeftCell = cells[leftColumn][row];
      topRightcell = cells[column][row];
      bottomRightCell = cells[column][bottomRow];
      bottomLeftCell = cells[leftColumn][bottomRow];
    }
    else {
      jAssert("getVelocityAt: Logic error in getting surrounding 4 cells.", false);
      topLeftCell = cells[column][row];
      topRightcell = cells[column][row];
      bottomRightCell = cells[column][row];
      bottomLeftCell = cells[column][row];
    }
    
    PVector topLeftVelocity = topLeftCell.getCenterVelocity();
    PVector topRightVelocity = topRightcell.getCenterVelocity();
    PVector bottomRightVelocity = bottomRightCell.getCenterVelocity();
    PVector bottomLeftVelocity = bottomLeftCell.getCenterVelocity();
  
    float velocityX = bilinearInterpolate(
      topLeftCell.position, topLeftVelocity.x, 
      topRightcell.position, topRightVelocity.x, 
      bottomRightCell.position, bottomRightVelocity.x, 
      bottomLeftCell.position, bottomLeftVelocity.x, 
      position);
      
    float velocityY = bilinearInterpolate(
      topLeftCell.position, topLeftVelocity.y, 
      topRightcell.position, topRightVelocity.y, 
      bottomRightCell.position, bottomRightVelocity.y, 
      bottomLeftCell.position, bottomLeftVelocity.y, 
      position);
      
    return new PVector(velocityX, velocityY);
  }


  /**
   * returns the column containing the given position
   * @param position
   *     MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  private int columnFromPosition(PVector position) {
//    jAssert("columnFromPosition: position.x is < 0: " + position, position.x >= 0);
//    jAssert("columnFromPosition: position.x is > regionWidth: " + position, position.x < cellWidth * cells.length);
    return constrain(floor(position.x / cellWidth), 0, cells.length-1);
  }

  /**
   * returns the row containing the given position
   * position MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  private int rowFromPosition(PVector position) {
//    jAssert("rowFromPosition: position.y is < 0: " + position, position.y >= 0);
//    jAssert("rowFromPosition: position.y is > regionWidth " + position, position.y < cellWidth * cells[0].length);
    return constrain(floor(position.y / cellWidth), 0, cells[0].length-1);
  }

  /**
   * Returns the cell containing the given position
   * position MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  private FluidGridCell getCellContaining(PVector position) {
    return cells[columnFromPosition(position)][rowFromPosition(position)];
  }
}

