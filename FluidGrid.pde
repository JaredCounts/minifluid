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
    assert(cellWidth > 0);
    assert(regionWidth > 0);
    assert(regionHeight > 0);

    this.cellWidth = cellWidth;

    int columnCount = ceil(regionWidth / cellWidth);
    int rowCount = ceil(regionHeight / cellWidth);

    cells = new FluidGridCell[columnCount][rowCount];

    for (int i = 0; i < columnCount; i++) {
      for (int j = 0; j < rowCount; j++) {
        cells[i][j] = new FluidGridCell();
        
        cells[i][j].position = new PVector(i * cellWidth + cellWidth/2.f, j * cellWidth + cellWidth/2.f);
        cells[i][j].topEdgePosition = new PVector(i * cellWidth + cellWidth / 2.f, j * cellWidth);
        cells[i][j].rightEdgePosition = new PVector(i * cellWidth + cellWidth, j * cellWidth + cellWidth/2.f);
        cells[i][j].bottomEdgePosition = new PVector(i * cellWidth + cellWidth / 2.f, j * cellWidth + cellWidth);
        cells[i][j].leftEdgePosition = new PVector(i * cellWidth, j * cellWidth + cellWidth/2.f);
        
        if (i != 0) { // initialize left/right edge velocities
          // both cells will point to the same velocity objects
          Float velocityX = new Float(0);
          cells[i][j].velocityXLeft = velocityX;
          cells[i-1][j].velocityXRight = velocityX;
        }
        
        if (j != 0) { // initialize top/bottom edge velocities
          Float velocityY = new Float(0);
          cells[i][j].velocityYTop = velocityY;
          cells[i][j-1].velocityYBottom = velocityY;
        }
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
    assert(timestep > 0);
     
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
    else
      solveCount = 1;
    
    for (int solve = 0; solve < solveCount; solve++) {
      /* -------------- External Forces -------------- */
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
      
      /* -------------- Convection -------------- */
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
      
      // to not contaminate velocities as we solve, store our new velocities into a temporary matrix
      // where newVelocities[i][j] gives us the velocity for upper/left wall of cell i,j
      // we size each by number of cells + 1 to account for very bottom/right edges
      float[][] newVelocitiesX = new float[cells.length+1][cells[0].length+1];
      float[][] newVelocitiesY = new float[cells.length+1][cells[0].length+1];
      for (int i = 0; i < cells.length; i++) {
        for (int j = 0; j < cells[i].length; j++) {
          FluidGridCell cell = cells[i][j];
          
          // solve for left edge
          PVector leftEdgePosition = cell.leftEdgePosition;
          PVector leftEdgeTracedPosition = traceFrom(leftEdgePosition, timestep); 
          PVector leftEdgeTracedVelocity = getVelocityAt(leftEdgeTracedPosition); 
          newVelocitiesX[i][j] = leftEdgeTracedVelocity.x;
          
          // solve for upper edge
          PVector topEdgePosition = cell.topEdgePosition;
          PVector topEdgeTracedPosition = traceFrom(topEdgePosition, timestep);
          PVector topEdgeTracedVelocity = getVelocityAt(topEdgeTracedPosition);
          newVelocitiesY[i][j] = topEdgeTracedVelocity.y;
          
          // lets not forget the last cells' bottom and right velocities
          if (i == cells.length-1) {
            PVector rightEdgePosition = cell.rightEdgePosition;
            PVector rightEdgeTracedPosition = traceFrom(rightEdgePosition, timestep);
            PVector rightEdgeTracedVelocity = getVelocityAt(rightEdgeTracedPosition);
            newVelocitiesX[i+1][j] = rightEdgeTracedVelocity.x;
          }
          if (j == cells[i].length-1) {
            PVector bottomEdgePosition = cell.bottomEdgePosition;
            PVector bottomEdgeTracedPosition = traceFrom(bottomEdgePosition, timestep);
            PVector bottomEdgeTracedVelocity = getVelocityAt(bottomEdgeTracedPosition);
            newVelocitiesY[i][j+1] = bottomEdgeTracedVelocity.y;
          }
        }
      }
      
      // now update our cell velocities
      for (int i = 0; i < cells.length; i++) {
        for (int j = 0; j < cells[i].length; j++) {
          FluidGridCell cell = cells[i][j];
          
          cell.velocityXLeft = newVelocitiesX[i][j];
          cell.velocityYTop  = newVelocitiesY[i][j];
          
          if (i == cells.length-1)
            cell.velocityXRight  = newVelocitiesX[i+1][j];
          if (j == cells[i].length-1)
            cell.velocityYBottom = newVelocitiesY[i][j+1];
        }
      }
      
      
      /* -------------- Viscosity -------------- */
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
            newVelocitiesX[i][j] = cell.velocityXLeft + viscousFactor * (-cell.velocityXLeft + cell.velocityXRight);
          }
          else {
            FluidGridCell leftCell = cells[i-1][j];
            newVelocitiesX[i][j] = cell.velocityXLeft + viscousFactor * (leftCell.velocityXLeft - 2f * cell.velocityXLeft + cell.velocityXRight);
          }
          
          // solve for upper edge
          if (j == 0) {
            newVelocitiesY[i][j] = cell.velocityYTop + viscousFactor * (-cell.velocityYTop + cell.velocityYBottom);
          }
          else {
            FluidGridCell topCell = cells[i][j-1];
            newVelocitiesY[i][j] = cell.velocityYTop + viscousFactor * (topCell.velocityYTop - 2f * cell.velocityYTop + cell.velocityYBottom);
          }
          
          // and lets not forget far right and far bottom edge velocities!
          if (i == cells.length-1) {
            newVelocitiesX[i+1][j] = cell.velocityXRight + viscousFactor * (cell.velocityXLeft - cell.velocityXRight);
          }
          if (i == cells.length-1) {
            newVelocitiesY[i][j+1] = cell.velocityYBottom + viscousFactor * (cell.velocityYTop - cell.velocityYBottom);
          }
          
        }
      }
      
      // now update our cell velocities
      // we should probably make this its own method
      for (int i = 0; i < cells.length; i++) {
        for (int j = 0; j < cells[i].length; j++) {
          FluidGridCell cell = cells[i][j];
          
          cell.velocityXLeft = newVelocitiesX[i][j];
          cell.velocityYTop  = newVelocitiesY[i][j];
          
          if (i == cells.length-1)
            cell.velocityXRight  = newVelocitiesX[i+1][j];
          if (j == cells[i].length-1)
            cell.velocityYBottom = newVelocitiesY[i][j+1];
        }
      }
      
      /* -------------- Incompressability and Pressure -------------- */
      // using successive over-relaxation
      //   Foster and Metaxas, "Realistic Animation of Liquids"
      // *NOTE* this is the easy approach. 
      // Fedkiw and Foster suggests using Preconditioning Conjugate Gradient (PCG) on the Laplacian linear system
      // relating incompressability and pressure
      // http://physbam.stanford.edu/~fedkiw/papers/stanford2001-02.pdf
      // idk what any of that means yet
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
    assert(finalPosition.x >= 0);
    assert(finalPosition.x <= cellWidth * cells.length);
    assert(finalPosition.y >= 0);
    assert(finalPosition.y <= cellWidth * cells[0].length);
    
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

    FluidGridCell topCell, bottomCell, leftCell, rightCell;
    if (topHalf) {
       topCell = cells[column][row-1];
       bottomCell = cells[column][row];
    }
    else {
       topCell = cells[column][row];
       bottomCell = cells[column][row+1];
    }
    if (leftHalf) {
      leftCell = cells[column-1][row];
      rightCell = cells[column][row];
    }
    else {
      leftCell = cells[column][row];
      rightCell = cells[column+1][row];
    }

    float velocityX = bilinearInterpolate(topCell.leftEdgePosition, topCell.velocityXLeft,
                                          topCell.rightEdgePosition, topCell.velocityXRight,
                                          bottomCell.rightEdgePosition, bottomCell.velocityXRight,
                                          bottomCell.leftEdgePosition, bottomCell.velocityXLeft,
                                          position);
    float velocityY = bilinearInterpolate(leftCell.topEdgePosition, leftCell.velocityYTop,
                                          leftCell.bottomEdgePosition, leftCell.velocityYBottom,
                                          rightCell.bottomEdgePosition, rightCell.velocityYBottom,
                                          rightCell.topEdgePosition, rightCell.velocityYTop,
                                          position);
                                          
    return new PVector(velocityX, velocityY);
  }


  /**
   * returns the column containing the given position
   * @param position
   *     MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  private int columnFromPosition(PVector position) {
    assert(position.x >= 0);
    assert(position.x <= cellWidth * cells.length);
    return floor(position.x / cellWidth);
  }

  /**
   * returns the row containing the given position
   * position MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  private int rowFromPosition(PVector position) {
    assert(position.y >= 0);
    assert(position.y <= cellWidth * cells[0].length);
    return floor(position.y / cellWidth);
  }

  /**
   * Returns the cell containing the given position
   * position MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  private FluidGridCell getCellContaining(PVector position) {
    return cells[columnFromPosition(position)][rowFromPosition(position)];
  }
}

