/** FluidGrid
 * MAC staggered grid
 * where velocities are stored at cell edges
 * and pressure in cell centers.
 * Mutable object.
 */
class FluidGrid {
  // cells will be square
  final float cellWidth;

  // grid of cells
  // to access, use cell cells[i][j], where cell is at the i'th column and j'th row.
  // cells[0][0]'s top left corner is at (0,0).
  private final FluidGridCell[][] cells;

  /** Constructor for fluid grid
   * Initializes all of the cells in the grid.
   * For best results, choose a cellWidth that's a common factor of regionWidth and regionHeight.
   * @param cellWidth height and width of each square cell
   *                  must be > 0
   * @param regionWidth width of the region to solve fluid in
   * @param regionHeight height of region
   */
  FluidGrid(float cellWidth, float regionWidth, float regionHeight) {
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
   * Solving
   * Using Foster and Fedkiw's "Practical Animation of Liquids" as reference
   * http://physbam.stanford.edu/~fedkiw/papers/stanford2001-02.pdf
   */
  public void solve() {
    // determine maximum timestep usinmg Courant-Friedrich's-Lewy (CFL) condition
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
    for (fluidGridCell[] cellColumn : cells) {
      for (fluidGridCell cell : cellColumn) {
        // at this point, incompressability should be solved for from the previous solve step
        // therefore the inward velocity of each cell should match the outward velocity
        // so all we need is one x value and one y value from the cell
        // another note, we don't need to do any expensive square-root operations until we found our max value
        float velocityMagnitudeSquared = sq(cell.velocityXLeft) + sq(cell.velocityYTop);
        maxVelocityMagitudeSquared = max(maxVelocityMagitudeSquared, velocityMagnitudeSquared);
      }
    }
    float maxVelocityMagnitude = sqrt(maxVelocityMagitudeSquared);
    
    // and finally, CFL
    float timestep = cellWidth / maxVelocityMagnitude;
     
    
    
    /* -------------- External Forces -------------- */
    // just integrate gravity into cell edge velocities
    for (int i = 0; i < fluidGridCell.length; i++) {
      for (int j = 0; j < fluidGridCell[i].length; j++) {
        FluidGridCell cell = fluidGridCell[i][j];
        // since edges are shared between cells, we only solve for left and top
        cell.velocityXLeft += GRAVITY.x * timestep;
        cell.velocityYTop += GRAVITY.y * timestep;
        
        // except for the bottom-most and right-most cells
        if (i == fluidGridCell.length-1)
          cell.velocityXRight += GRAVITY.x * timestep;
        if (j == fluidGridCell[i].length-1)
          cell.velocityYBottom += GRAVITY.y * timestep;
      }
    }
    
    /* -------------- Convection -------------- */
    // using semi-lagrangian method
    //   Stam, "Stable Fluids"
    //   http://www.autodeskresearch.com/pdf/ns.pdf
    
    
    /* -------------- Viscosity -------------- */
    // using standard central differencing
    //   Foster and Metaxas, "Realistic Animation of Liquids"
    //   http://graphics.stanford.edu/courses/cs468-05-fall/Papers/foster-metaxas-gmip96.pdf
    
    
    /* -------------- Incompressability and Pressure -------------- */
    // using successive over-relaxation
    //   Foster and Metaxas, "Realistic Animation of Liquids"
    // *NOTE* this is the easy approach. 
    // Fedkiw and Foster suggests using Preconditioning Conjugate Gradient (PCG) on the Laplacian linear system
    // relating incompressability and pressure
    // http://physbam.stanford.edu/~fedkiw/papers/stanford2001-02.pdf
    // idk what any of that means yet

  }
  
  /**
   * Returns a velocity at the given position
   */
  public PVector getVelocityAt(PVector position) {
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
   * returns the column containing position
   * position MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  int columnFromPosition(PVector position) {
    assert(position.x >= 0);
    assert(position.x <= cellWidth * cells.length);
    return floor(position.x / cellWidth);
  }

  /**
   * returns the row containing position
   * position MUST BE within (0,0) -> (regionWidth, regionHeight) rectangle
   */
  int rowFromPosition(PVector position) {
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

