/** 
 * Fluid Grid Cell
 * Mutable object which acts as a container for pressure, and has references to its edge velocities
 */
 public enum CellType {
  SOLID,
  AIR,
  LIQUID
}
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
  
  // determine whether this cell is SOLID, AIR, or LIQUID
  CellType type;

  /**
   * Construct an empty FluidGridCell
   * note that edge velocities must be set manually. 
   */
  FluidGridCell(CellType cellType) {
    pressure = 0;
    type = cellType;
  }
  
  void draw() {
    PVector velocity = getCenterVelocity();
    colorMode(HSB,255);
    if (isLiquid())
      fill(pressure*100, 255, 255);
    else if (isSolid())
      fill(0);
    else
      fill(255);
    
    noStroke();
    rect(position.x-cellWidth/2, position.y-cellWidth/2, cellWidth, cellWidth);
    
    stroke(0);
    line(position.x, position.y, // from
         position.x + velocity.x, position.y + velocity.y); // to
  }
  
  PVector getCenterVelocity() {
    return new PVector((velocityXLeft + velocityXRight) / 2.f,
                                   (velocityYTop + velocityYBottom) / 2.f);
  }
  
  boolean isLiquid() {
    return type == CellType.LIQUID; 
  }
  
  boolean isSolid() {
    return type == CellType.SOLID; 
  }
  
  boolean isAir() {
    return type == CellType.AIR;
  }
  
}