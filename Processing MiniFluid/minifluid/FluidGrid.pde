/** FluidGrid
 * MAC staggered grid
 * where velocities are stored at cell edges
 * and pressure in cell centers.
 * Mutable object.
 */
class FluidGrid {
  // cells will be square
  final double cellWidth;
  
  // grid of cells
  // to access, use cell cells[i][j], where cell is at the i'th column and j'th row.
  final FluidGridCell[][] cells;
  
  /** Constructor for fluid grid
   * Initializes all of the cells in the grid.
   * For best results, choose a cellWidth that's a common factor of regionWidth and regionHeight.
   * @param cellWidth height and width of each square cell
   *                  must be > 0
   * @param regionWidth width of the region to solve fluid in
   * @param regionHeight height of region
   */
  FluidGrid(double cellWidth, double regionWidth, double regionHeight) {
    assert(cellWidth > 0);
    assert(regionWidth > 0);
    assert(regionHeight > 0);
    
    this.cellWidth = cellWidth;
    
    int columnCount = (int)Math.ceil(regionWidth / cellWidth);
    int rowCount = (int)Math.ceil(regionHeight / cellWidth);
    
    cells = new FluidGridCell[columnCount][rowCount];
    
    for (int i = 0; i < columnCount; i++) {
      for (int j = 0; j < rowCount; j++) {
        cells[i][j] = new FluidGridCell();
        if (i != 0) { // initialize left/right edge velocities
          // both cells will point to the same velocity objects
          Double velocityX = new Double(0);
          cells[i][j].velocityXLeft = velocityX;
          cells[i-1][j].velocityXRight = velocityX; 
        }
        if (j != 0) { // initialize top/bottom edge velocities
          Double velocityY = new Double(0);
          cells[i][j].velocityYTop = velocityY;
          cells[i][j-1].velocityYBottom = velocityY;
        }
      } 
    }
  }
}
/** FluidGridCell
 * Mutable object which acts as a container for pressure, and has references to its edge velocities
 */
class FluidGridCell {
  // pressure at the center of this cell
  double pressure; 
  
  // velocity objects which are shared with adjacent cells
  // must set these manually
  Double velocityXLeft, velocityXRight, velocityYTop, velocityYBottom;
  
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
