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
/** FluidGridCell
 * Mutable object which acts as a container for pressure, and has references to its edge velocities
 */
class FluidGridCell {
  // pressure at the center of this cell
  float pressure; 

  // velocity objects which are shared with adjacent cells
  // must set these manually
  Float velocityXLeft, velocityXRight, velocityYTop, velocityYBottom;

  // must set positions manually
  // position is the center of the cell
  PVector position;
  // each edge position is the center of each edge of the cell
  PVector topEdgePosition, rightEdgePosition, bottomEdgePosition, leftEdgePosition;
  
  // flag for whether this cell is filled with liquid
  boolean hasLiquid;

  /**
   * Construct an empty FluidGridCell
   * note that edge velocities must be set manually. 
   */
  FluidGridCell() {
    pressure = 0;
    hasLiquid = false;
  }
 
}

